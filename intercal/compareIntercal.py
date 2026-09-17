#!/usr/bin/env python3
# ---------------------------------------------------------------------------
#  Compare the MIP peak bar by bar across three intercalibration levels.
#  Produces three summary plots (L, R, L-R), a text table, and - unless
#  disabled - one control plot per individual fit.
#  Depends on ROOT only.
#
#
#  For running: unset DISPLAY
#               python3 macros/compareIntercal.py
#
#  Usage:
#    python3 macros/compareIntercal.py
#    python3 macros/compareIntercal.py --mode landau --fracHi 0.20
#    python3 macros/compareIntercal.py --noFitPlots
# ---------------------------------------------------------------------------
import os, argparse
import ROOT
 
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
 
p = argparse.ArgumentParser(description="MIP peak vs bar for three intercalibration levels")
p.add_argument("--plotDir",    default="/eos/user/t/tzini/TestBeamAnalysis/plot")
p.add_argument("--inputDir",   default="/afs/cern.ch/user/t/tzini/private/tzini/Lab5015Analysis/plots")
p.add_argument("--vov",        type=float, default=3.00, help="over-voltage")
p.add_argument("--th",         type=int,   default=10,   help="threshold")
p.add_argument("--searchMin",  type=float, default=350., help="minimum energy where the peak search starts")
p.add_argument("--nSmooth",    type=int,   default=4,    help="half-width of the moving average, in bins")
p.add_argument("--relLo",      type=float, default=0.30, help="left edge of the fit window, as a fraction of the peak position")
p.add_argument("--relHi",      type=float, default=0.35, help="right edge of the fit window, as a fraction of the peak position")
p.add_argument("--noFitPlots", action="store_true",      help="skip the per-fit control plots")
a = p.parse_args()
 
# --- the three configurations: label, short tag, file, colour, marker ------
#     colours taken from the CMS palette already used in macros/utils.py
CONFIGS = [
    ("No intercalibration",     "NoCal",     "moduleCharacterization_step1_3880_NoCal.root",     "#3f90da", 20),
    ("Legacy intercalibration", "LegacyCal", "moduleCharacterization_step1_3880_LegacyCal.root", "#e76300", 21),
    ("New intercalibration",    "NewCal",    "moduleCharacterization_step1_NewCal.root",         "#bd1f01", 22),
]
SIDES = ["L", "R", "L-R"]
NBARS = 16
 
outdir = os.path.join(a.plotDir, "Intercal_plots")
os.makedirs(outdir, exist_ok=True)
 
 
def fitPeak(h, savePath=None, titleLine=""):
    """Find the MIP peak and fit it with a Landau.
       Returns the maximum of the fitted curve, or None when unreliable.
       If savePath is given, writes a control plot of this single fit."""
    if not h or h.GetEntries() < 200:
        return None
    nb = h.GetNbinsX()
    nSm = a.nSmooth
    b1 = max(h.FindBin(a.searchMin), 1 + nSm)
 
    # --- locate the peak: maximum of the moving average, which is stable even
    #     on a flat top, unlike the single highest bin
    best, bpk = -1., b1
    for b in range(b1, nb - nSm + 1):
        s = sum(h.GetBinContent(b + k) for k in range(-nSm, nSm + 1))
        if s > best:
            best, bpk = s, b
    if best <= 0.:
        return None
 
    ymax = best / (2. * nSm + 1.)
    xpk = h.GetBinCenter(bpk)
 
    # --- fit window: a fixed span around the peak, wide enough to contain the
    #     whole shape including the right-hand tail
    blo = max(h.FindBin((1. - a.relLo) * xpk), 1)
    bhi = min(h.FindBin((1. + a.relHi) * xpk), nb)
    if bhi - blo < 5:
        return xpk
 
    xlo, xhi = h.GetBinCenter(blo), h.GetBinCenter(bhi)
 
    # --- Landau fit
    f = ROOT.TF1("f_" + h.GetName(), "[0]*TMath::Landau(x,[1],[2])", xlo, xhi)
    f.SetParameters(ymax, xpk, 0.10 * xpk)
    f.SetParLimits(1, xlo, xhi)
    f.SetParLimits(2, 0.01 * xpk, 0.50 * xpk)
    h.Fit(f, "QRN")
 
    # the Landau maximum sits at mpv - 0.22278*sigma, so the peak is not
    # GetParameter(1): take the maximum of the curve itself
    xpeak = f.GetMaximumX(xlo, xhi)
 
    ok = (xlo < xpeak < xhi)
    if not ok:
        xpeak = xpk                                  # fall back on the bin maximum
 
    # --- control plot of this single fit ----------------------------------
    if savePath:
        c = ROOT.TCanvas("cfit", "cfit", 800, 600)
        c.SetLeftMargin(0.13)
        c.SetBottomMargin(0.13)
        h.SetLineColor(ROOT.kGray + 2)
        h.SetLineWidth(1)
        h.GetXaxis().SetTitle("energy [ADC]")
        h.GetYaxis().SetTitle("entries")
        h.GetXaxis().SetTitleSize(0.05); h.GetXaxis().SetLabelSize(0.045)
        h.GetYaxis().SetTitleSize(0.05); h.GetYaxis().SetLabelSize(0.045)
        h.GetYaxis().SetTitleOffset(1.25)
        h.Draw("hist")
 
        f.SetLineColor(ROOT.TColor.GetColor("#bd1f01"))
        f.SetLineWidth(2)
        f.Draw("same")
 
        ln = ROOT.TLine(xpeak, 0., xpeak, 1.05 * h.GetMaximum())
        ln.SetLineColor(ROOT.TColor.GetColor("#3f90da"))
        ln.SetLineStyle(2); ln.SetLineWidth(2)
        ln.Draw("same")
 
        t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.038)
        t.DrawLatex(0.52, 0.86, titleLine)
        t.DrawLatex(0.52, 0.80, "peak = %.1f ADC" % xpeak)
        t.DrawLatex(0.52, 0.74, "fit range [%.0f, %.0f]" % (xlo, xhi))
        if not ok:
            t.SetTextColor(ROOT.TColor.GetColor("#bd1f01"))
            t.DrawLatex(0.52, 0.68, "fit outside range - bin used")
 
        c.Print(savePath)
        del c
 
    return xpeak
 
 
# --- peak extraction -------------------------------------------------------
peaks = {}          # (label, side, bar) -> peak energy
for label, tag, fname, _, _ in CONFIGS:
    path = os.path.join(a.inputDir, fname)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print("[ERROR] cannot open %s" % path)
        continue
 
    fitdir = None
    if not a.noFitPlots:
        fitdir = os.path.join(outdir, "fit_check", tag)
        os.makedirs(fitdir, exist_ok=True)
 
    for side in SIDES:
        for b in range(NBARS):
            hname = "h1_energy_bar%02d%s_Vov%.2f_th%02d" % (b, side, a.vov, a.th)
            h = f.Get(hname)
            if not h:
                continue
            h = h.Clone(hname + "_c")
            h.SetDirectory(0)
 
            savePath = None
            if fitdir:
                savePath = "%s/fit_bar%02d%s.png" % (fitdir, b, side.replace("-", ""))
            v = fitPeak(h, savePath, "%s - bar%02d %s" % (label, b, side))
            if v:
                peaks[(label, side, b)] = v
    f.Close()
    print("[INFO] read %s" % fname)
 
if not peaks:
    raise SystemExit("[ERROR] no peak extracted: check the file names and --vov / --th")
 
 
# --- relative spread, the number that summarises each curve ----------------
def spread(label, side):
    v = [peaks[(label, side, b)] for b in range(NBARS) if (label, side, b) in peaks]
    if len(v) < 2:
        return None, None
    m = sum(v) / len(v)
    rms = (sum((x - m) ** 2 for x in v) / len(v)) ** 0.5
    return m, 100. * rms / m
 
 
# --- one summary plot per side --------------------------------------------
for side in SIDES:
    c = ROOT.TCanvas("c_" + side, "c_" + side, 900, 650)
    c.SetGridy()
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.13)
 
    mg = ROOT.TMultiGraph()
    leg = ROOT.TLegend(0.16, 0.74, 0.64, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.033)
 
    keep = []
    for label, _, _, hexcol, marker in CONFIGS:
        col = ROOT.TColor.GetColor(hexcol)
        g = ROOT.TGraph()
        for b in range(NBARS):
            if (label, side, b) in peaks:
                g.SetPoint(g.GetN(), b, peaks[(label, side, b)])
        if g.GetN() == 0:
            continue
        g.SetLineColor(col);   g.SetLineWidth(2)
        g.SetMarkerColor(col); g.SetMarkerStyle(marker); g.SetMarkerSize(1.3)
        mg.Add(g, "PL")
        keep.append(g)
 
        m, rel = spread(label, side)
        leg.AddEntry(g, "%s   #LTE#GT = %.0f,  RMS = %.1f%%" % (label, m, rel), "lp")
 
    mg.Draw("A")
    mg.GetXaxis().SetTitle("bar number")
    mg.GetYaxis().SetTitle("MIP peak [ADC]")
    mg.GetXaxis().SetTitleSize(0.05); mg.GetXaxis().SetLabelSize(0.045)
    mg.GetYaxis().SetTitleSize(0.05); mg.GetYaxis().SetLabelSize(0.045)
    mg.GetYaxis().SetTitleOffset(1.25)
    mg.GetXaxis().SetLimits(-0.5, NBARS - 0.5)
    leg.Draw()
 
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.032); t.SetTextFont(42)
    t.DrawLatex(0.63, 0.92, "side %s   V_{ov} = %.2f   threshold %d" % (side, a.vov, a.th))
 
    name = side.replace("-", "")
    for ext in ("png", "pdf"):
        c.Print("%s/mipPeak_vs_bar_%s.%s" % (outdir, name, ext))
    del c
 
# --- the same information as a table --------------------------------------
with open(os.path.join(outdir, "mipPeak_vs_bar.txt"), "w") as out:
    out.write("# MIP peak [ADC] from a Landau fit, Vov %.2f, threshold %d, window [-%.0f%%, +%.0f%%]\n"
              % (a.vov, a.th, 100. * a.relLo, 100. * a.relHi))
    for side in SIDES:
        out.write("\n# side %s\n" % side)
        out.write("%-5s" % "bar" + "".join("%26s" % l for l, _, _, _, _ in CONFIGS) + "\n")
        for b in range(NBARS):
            row = "%-5d" % b
            for label, _, _, _, _ in CONFIGS:
                v = peaks.get((label, side, b))
                row += "%26s" % ("%.1f" % v if v else "-")
            out.write(row + "\n")
        for label, _, _, _, _ in CONFIGS:
            m, rel = spread(label, side)
            if m:
                out.write("# %-26s mean %.1f   RMS %.1f%%\n" % (label, m, rel))
 
print("\n[INFO] plots and table written to: %s\n" % outdir)
for side in SIDES:
    for label, _, _, _, _ in CONFIGS:
        m, rel = spread(label, side)
        if m:
            print("  %-5s %-26s mean %7.1f   RMS %5.1f%%" % (side, label, m, rel))
 
#!/usr/bin/env python3
# ---------------------------------------------------------------------------
#  Correlation between two sets of intercalibration factors, channel by channel.
#
#  Reads two csv files written by writeTOFHIRcalib.py (columns
#  bar,side,vov,th,calib), matches them on (bar, side) at one given (vov, th),
#  and produces:
#    - a scatter plot, one point per channel, with the y = x reference and a
#      straight-line fit;
#    - a residual plot (new - reference) channel by channel, where structure
#      is far easier to see than in the scatter;
#    - a text summary with correlation, slope, spreads and the worst channels.
#
#  Depends on ROOT only.
#
#  Usage, from the repository root:
#    python3 intercal/plotCalibCorrelation.py \
#        --refFile intercal/CONF76/TOFHIR_LO_calibration_factors.csv \
#        --newFile intercal/CONF3880/TOFHIR_LO_calibration_factors.csv
# ---------------------------------------------------------------------------
import os, csv, argparse
import ROOT
 
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
 
p = argparse.ArgumentParser(description="Correlation between two sets of intercalibration factors")
p.add_argument("--refFile",  required=True, help="csv taken as reference (x axis)")
p.add_argument("--newFile",  required=True, help="csv to be compared (y axis)")
p.add_argument("--refLabel", default="Legacy", help="axis label for the reference set")
p.add_argument("--newLabel", default="New",    help="axis label for the compared set")
p.add_argument("--vov",      type=float, default=3.00)
p.add_argument("--th",       type=int,   default=10)
p.add_argument("--outDir",   default="/eos/user/t/tzini/TestBeamAnalysis/plot/Intercal_plots")
p.add_argument("--nLabel",   type=int, default=4, help="how many outlying channels to label directly")
a = p.parse_args()
 
# side -> colour, marker. Colours from the CMS palette used in macros/utils.py;
# the marker differs too, so identity never rests on colour alone.
STYLE = {"L": ("#3f90da", 20), "R": ("#e76300", 21)}
 
os.makedirs(a.outDir, exist_ok=True)
 
 
def load(path):
    """(bar, side) -> calib, for the requested (vov, th)."""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if abs(float(r["vov"]) - a.vov) > 1e-6 or int(r["th"]) != a.th:
                continue
            out[(int(r["bar"]), r["side"])] = float(r["calib"])
    return out
 
 
ref, new = load(a.refFile), load(a.newFile)
keys = sorted(set(ref) & set(new))
if not keys:
    raise SystemExit("[ERROR] no channel in common at Vov %.2f, threshold %d" % (a.vov, a.th))
 
only_ref = sorted(set(ref) - set(new))
only_new = sorted(set(new) - set(ref))
if only_ref or only_new:
    print("[WARNING] %d channels only in the reference, %d only in the new set"
          % (len(only_ref), len(only_new)))
 
x = [ref[k] for k in keys]
y = [new[k] for k in keys]
n = len(keys)
 
 
def mean(v):
    return sum(v) / len(v)
 
 
def rms(v):
    m = mean(v)
    return (sum((t - m) ** 2 for t in v) / len(v)) ** 0.5
 
 
def linfit(u, v):
    """Ordinary least squares of v on u. Returns (slope, intercept, pearson).
       Note this is not symmetric: it assumes u carries no error. With both
       axes measured the slope is biased low, so check that the bias cannot
       account for a slope far from 1 before reading anything into it."""
    if len(u) < 2:
        return 0., 0., 0.
    mu, mv = mean(u), mean(v)
    suv = sum((i - mu) * (j - mv) for i, j in zip(u, v))
    suu = sum((i - mu) ** 2 for i in u)
    svv = sum((j - mv) ** 2 for j in v)
    sl = suv / suu if suu > 0 else 0.
    return sl, mv - sl * mu, (suv / (suu * svv) ** 0.5 if suu * svv > 0 else 0.)
 
 
slope, inter, pearson = linfit(x, y)
 
# the two sides behave differently, so each also gets its own fit
sideFit = {}
for side in STYLE:
    xs = [ref[k] for k in keys if k[1] == side]
    ys = [new[k] for k in keys if k[1] == side]
    if len(xs) >= 2:
        sideFit[side] = linfit(xs, ys)
 
resid = [v - u for u, v in zip(x, y)]
 
# ---------------------------------------------------------------- scatter --
c = ROOT.TCanvas("c_corr", "c_corr", 750, 750)
c.SetLeftMargin(0.14); c.SetRightMargin(0.05)
c.SetBottomMargin(0.13); c.SetTopMargin(0.06)
c.SetGridx(); c.SetGridy()
 
lo = min(min(x), min(y)); hi = max(max(x), max(y))
pad = 0.08 * (hi - lo)
lo -= pad; hi += pad
 
frame = ROOT.TH2F("frame", "", 10, lo, hi, 10, lo, hi)
frame.GetXaxis().SetTitle("%s calibration factor" % a.refLabel)
frame.GetYaxis().SetTitle("%s calibration factor" % a.newLabel)
frame.GetXaxis().SetTitleSize(0.045); frame.GetXaxis().SetLabelSize(0.04)
frame.GetYaxis().SetTitleSize(0.045); frame.GetYaxis().SetLabelSize(0.04)
frame.GetYaxis().SetTitleOffset(1.45)
frame.Draw()
 
# y = x reference, in neutral ink: it is context, not a series
diag = ROOT.TLine(lo, lo, hi, hi)
diag.SetLineColor(ROOT.kGray + 1); diag.SetLineStyle(2); diag.SetLineWidth(2)
diag.Draw("same")
 
# one straight-line fit per side, in the colour of its own points. The fit
# over all channels together is quoted as a number instead of a fourth line:
# four lines in one square make the figure harder to read, not richer.
fitLines = {}
for side, (hexcol, _) in STYLE.items():
    if side not in sideFit:
        continue
    sl, it, _ = sideFit[side]
    ln = ROOT.TLine(lo, sl * lo + it, hi, sl * hi + it)
    ln.SetLineColor(ROOT.TColor.GetColor(hexcol))
    ln.SetLineWidth(2)
    ln.Draw("same")
    fitLines[side] = ln
 
graphs = {}
for side, (hexcol, marker) in STYLE.items():
    g = ROOT.TGraph()
    for k in keys:
        if k[1] == side:
            g.SetPoint(g.GetN(), ref[k], new[k])
    if g.GetN() == 0:
        continue
    col = ROOT.TColor.GetColor(hexcol)
    g.SetMarkerColor(col); g.SetMarkerStyle(marker); g.SetMarkerSize(1.4)
    g.SetLineColor(col)
    g.Draw("P same")
    graphs[side] = g
 
leg = ROOT.TLegend(0.17, 0.75, 0.40, 0.88)
leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.033)
for side, g in graphs.items():
    leg.AddEntry(g, "side %s" % side, "lp")
leg.AddEntry(diag, "y = x", "l")
leg.Draw()
 
# label only the few channels furthest from y = x
order = sorted(range(n), key=lambda i: -abs(resid[i]))
lab = ROOT.TLatex(); lab.SetTextFont(42); lab.SetTextSize(0.028)
lab.SetTextColor(ROOT.kGray + 3)
for i in order[:a.nLabel]:
    lab.DrawLatex(x[i] + 0.012 * (hi - lo), y[i], "bar%02d%s" % keys[i])
 
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.030)
t.DrawLatex(0.48, 0.30, "n = %d channels" % n)
ypos = 0.26
for side in sorted(sideFit):
    sl, _, rr = sideFit[side]
    t.DrawLatex(0.48, ypos, "side %s:  slope = %.3f,  #rho = %.3f" % (side, sl, rr))
    ypos -= 0.04
t.DrawLatex(0.48, ypos, "RMS of (%s - %s) = %.3f" % (a.newLabel, a.refLabel, rms(resid)))
t.SetTextSize(0.028)
t.DrawLatex(0.50, 0.93, "V_{ov} = %.2f   threshold %d" % (a.vov, a.th))
 
for ext in ("png", "pdf"):
    c.Print("%s/calib_correlation.%s" % (a.outDir, ext))
del c
 
# -------------------------------------------------------------- residuals --
c2 = ROOT.TCanvas("c_res", "c_res", 1000, 600)
c2.SetLeftMargin(0.12); c2.SetRightMargin(0.04)
c2.SetBottomMargin(0.16); c2.SetTopMargin(0.06)
c2.SetGridy()
 
amp = max(abs(min(resid)), abs(max(resid))) * 1.35
f2 = ROOT.TH1F("f2", "", n, -0.5, n - 0.5)
f2.GetYaxis().SetRangeUser(-amp, amp)
f2.GetYaxis().SetTitle("%s - %s" % (a.newLabel, a.refLabel))
f2.GetYaxis().SetTitleSize(0.05); f2.GetYaxis().SetLabelSize(0.045)
f2.GetYaxis().SetTitleOffset(1.05)
f2.GetXaxis().SetLabelSize(0.045)
for i, k in enumerate(keys):
    f2.GetXaxis().SetBinLabel(i + 1, "%d%s" % k)
f2.GetXaxis().LabelsOption("v")
f2.Draw()
 
zero = ROOT.TLine(-0.5, 0., n - 0.5, 0.)
zero.SetLineColor(ROOT.kGray + 2); zero.SetLineWidth(2)
zero.Draw("same")
 
gres = {}
for side, (hexcol, marker) in STYLE.items():
    g = ROOT.TGraph()
    for i, k in enumerate(keys):
        if k[1] == side:
            g.SetPoint(g.GetN(), i, resid[i])
    if g.GetN() == 0:
        continue
    col = ROOT.TColor.GetColor(hexcol)
    g.SetMarkerColor(col); g.SetMarkerStyle(marker); g.SetMarkerSize(1.4)
    g.Draw("P same")
    gres[side] = g
 
leg2 = ROOT.TLegend(0.14, 0.80, 0.34, 0.91)
leg2.SetBorderSize(0); leg2.SetFillStyle(0); leg2.SetTextSize(0.040)
for side, g in gres.items():
    leg2.AddEntry(g, "side %s" % side, "p")
leg2.Draw()
 
t2 = ROOT.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.038)
t2.DrawLatex(0.60, 0.88, "mean %+.4f,  RMS %.4f" % (mean(resid), rms(resid)))
t2.SetTextSize(0.032)
t2.DrawLatex(0.60, 0.93, "V_{ov} = %.2f   threshold %d" % (a.vov, a.th))
 
for ext in ("png", "pdf"):
    c2.Print("%s/calib_residuals.%s" % (a.outDir, ext))
del c2
 
# ------------------------------------------------------------ text output --
txt = os.path.join(a.outDir, "calib_correlation.txt")
with open(txt, "w") as out:
    out.write("# %s (x) vs %s (y), Vov %.2f, threshold %d\n" % (a.refLabel, a.newLabel, a.vov, a.th))
    out.write("# reference : %s\n# compared  : %s\n\n" % (a.refFile, a.newFile))
    out.write("channels matched      %d\n" % n)
    out.write("Pearson correlation   %.4f\n" % pearson)
    out.write("slope / intercept     %.4f / %.4f   (all channels)\n" % (slope, inter))
    for side in sorted(sideFit):
        sl, it, rr = sideFit[side]
        out.write("slope / intercept     %.4f / %.4f   (side %s, rho %.4f)\n" % (sl, it, side, rr))
    out.write("spread (RMS) of %-8s %.4f\n" % (a.refLabel, rms(x)))
    out.write("spread (RMS) of %-8s %.4f\n" % (a.newLabel, rms(y)))
    out.write("difference: mean %+.4f   RMS %.4f   min %+.4f   max %+.4f\n\n"
              % (mean(resid), rms(resid), min(resid), max(resid)))
    out.write("%-8s %10s %10s %10s\n" % ("channel", a.refLabel, a.newLabel, "diff"))
    for i in order:
        out.write("bar%02d%-3s %10.4f %10.4f %+10.4f\n" % (keys[i][0], keys[i][1], x[i], y[i], resid[i]))
 
print("\n  channels matched    %d" % n)
print("  Pearson             %.4f" % pearson)
print("  slope / intercept   %.4f / %.4f   (all channels)" % (slope, inter))
for side in sorted(sideFit):
    sl, it, rr = sideFit[side]
    print("  slope / intercept   %.4f / %.4f   (side %s, rho %.4f)" % (sl, it, side, rr))
print("  RMS %-8s        %.4f" % (a.refLabel, rms(x)))
print("  RMS %-8s        %.4f" % (a.newLabel, rms(y)))
print("  difference          mean %+.4f   RMS %.4f" % (mean(resid), rms(resid)))
print("\n  worst channels:")
for i in order[:a.nLabel]:
    print("    bar%02d%-3s %s %.3f   %s %.3f   diff %+.3f"
          % (keys[i][0], keys[i][1], a.refLabel, x[i], a.newLabel, y[i], resid[i]))
print("\n[INFO] plots and summary written to: %s\n" % a.outDir)
 
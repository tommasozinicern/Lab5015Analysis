void giveSens(std::string fileName = "/afs/cern.ch/user/t/tzini/private/tzini/Lab5015Analysis/plots/moduleCharacterization_step2_generalLabel.root",
             float Vov = 3.00, int vth = 10)  
             
             // inserire qui il path ai dati dello step2, l OV e la vth
{





  gStyle->SetOptFit(1111); // crea la tabella con i vari valori del fit
  gStyle->SetOptStat(0);   //taglia il riquadro che root fa di d efault

  TFile* f = TFile::Open(fileName.c_str());  //apertura del file e controllo
  if(!f || f->IsZombie()){
    std::cout << "Impossibile aprire il file " << fileName << std::endl;
    return;
  }

  // pattern del nome degli istogrammi finali (energyRatio+totRatio+fase corrette)
  TString prefix = "h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_bar";
  TString tag = Form("_Vov%.2f_th%02d_", Vov, vth);

  TIter next(f->GetListOfKeys());  //preparazione dello scorrimento degli oggetti nel file
  TKey* key;

  std::map<int, std::pair<double,double>> results; // bar -> (risoluzione barra [ps], errore)

  TCanvas* c = new TCanvas("c_CTR", "c_CTR", 800, 600);

  while( (key = (TKey*)next()) )  //inizio ciclo
  
  {
    TString name = key->GetName();               //se trovo il file corretto si salta subito alla prox iterazione del ciclo
    if( !name.BeginsWith(prefix) ) continue;
    if( !name.Contains(tag) ) continue;

    TH1F* h = (TH1F*)f->Get(name);
    if( !h || h->GetEntries() < 50 ) continue;   // salta istogrammi troppo vuoti/vuoti

    // estrae il numero di barra dai due caratteri dopo "bar" nel nome
    int barPos = name.Index("bar") + 3;
    TString sub = name(barPos, 2);
    int barID = sub.Atoi();

    if( results.find(barID) != results.end() ) continue; // già trovato (primo energyBin disponibile)

    // fit gaussiano intorno al picco della distribuzione
    double mean = h->GetBinCenter(h->GetMaximumBin());
    double rms  = h->GetRMS();
    TF1* fitFunc = new TF1("fitFunc", "gaus", mean - 1.5*rms, mean + 1.5*rms); //provare a 1rms (o adattare a metà altezza gaus) -> trovare valore int.......................................................
    h->Fit(fitFunc, "QR");

    double sigma    = fitFunc->GetParameter(2);  //estraggo i dati relativi alla sensibilità
    double sigmaErr = fitFunc->GetParError(2);

    // la risoluzione della singola barra e' meta' della sigma di deltaT = tR - tL
    double barRes    = sigma / 2.;
    double barResErr = sigmaErr / 2.;
    results[barID] = std::make_pair(barRes, barResErr);

    h->SetTitle(Form("bar %02d, Vov=%.2f, th=%d;#Delta t = t_{R}-t_{L}  [ps];eventi", barID, Vov, vth));
    h->Draw();
    c->Print(Form("CTR_fit_bar%02d.png", barID));
  }

  std::cout << "\n=== Risoluzione temporale per barra (Vov=" << Vov << ", th=" << vth << ") ===\n";
  std::cout << "bar\trisoluzione barra [ps]\n";

  TGraphErrors* gr = new TGraphErrors();  //tabella e grafico riassuntivo
  int ip = 0;
  for(auto& r : results)
  {
    std::cout << r.first << "\t" << r.second.first << " +- " << r.second.second << std::endl;
    gr->SetPoint(ip, r.first, r.second.first);
    gr->SetPointError(ip, 0, r.second.second);
    ip++;
  }

  TCanvas* c2 = new TCanvas("c_summary", "c_summary", 800, 600);
  gr->SetTitle(Form("risoluzione temporale vs barra (Vov=%.2f, th=%d);bar ID;risoluzione [ps]", Vov, vth));
  gr->SetMarkerStyle(20);
  gr->Draw("AP");
  gr->GetYaxis()->SetRangeUser(0, 40);
  c2->Print("CTR_summary_vs_bar.png");

  std::cout << "\nSalvati: CTR_fit_barXX.png per ogni barra, e CTR_summary_vs_bar.png col riepilogo.\n";
}
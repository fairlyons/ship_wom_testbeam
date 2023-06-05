void cell_ana(TString f) {
    gStyle->SetOptStat("eou");

    TChain *chain = new TChain("Detected");
    chain->Add(f);

    TTreeReader *data = new TTreeReader(chain);
    TTreeReaderValue<Int_t>      wom(*data,"sipmNumber");
    TTreeReaderValue<Double_t>   time(*data,"time");
    TTreeReaderValue<Int_t>      eventN(*data,"eventNumber");

    int entries = chain->GetEntries();
    int events = 1;

    TH1I *h1 = new TH1I("", "", 16, 0, 16);
    h1->GetXaxis()->SetTitle("SiPM Number");
    h1->GetYaxis()->SetTitle("Events");

    TH1D *h2 = new TH1D("", "", 100, 0, 1000);
    h2->GetXaxis()->SetTitle("Time (ns)");
    h2->GetYaxis()->SetTitle("Events");

    double N[events], up[events], down[events];
    for(int i = 0; i < events; i++) N[i] = i+1;

    while (data->Next()) {
      h1->Fill(*wom);
      h2->Fill(*time);
      if(*wom < 8) down[*eventN]++;
      if(*wom > 7) up[*eventN]++;
    }

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,800);
    //gPad->SetLogy();
    h1->SetLineColor(kBlue);
    h1->Draw("HIST");

    TCanvas *c2 = new TCanvas("c2","c2",10,10,1200,800);
    //gPad->SetLogy();
    h2->SetLineColor(kBlue);
    h2->Draw("HIST");

    TCanvas *c3 = new TCanvas("c3","c3",10,10,1200,800);
    c3->Divide(1,2);
    c3->cd(1);
    TGraph *h3 = new TGraph(events, N, up);
    h3->GetXaxis()->SetTitle("Event");
    h3->GetYaxis()->SetTitle("Detected photons in Upper WOM");
    //gPad->SetLogy();
    h3->SetFillColor(38);
    h3->Draw("AB");
    c3->cd(2);
    TGraph *h4 = new TGraph(events, N, down);
    h4->GetXaxis()->SetTitle("Event");
    h4->GetYaxis()->SetTitle("Detected photons in Lower WOM");
    //gPad->SetLogy();
    h4->SetFillColor(38);
    h4->Draw("AB");

    return;
}

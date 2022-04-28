void cell_ana() {
    gStyle->SetOptStat("eou");

    TString f = "b/cerenkov.root";
    TChain *chain = new TChain("Detected");
    chain->Add(f);

    TTreeReader *data = new TTreeReader(chain);
    TTreeReaderValue<Double_t>   x(*data,"row_wise_branch.x");
    TTreeReaderValue<Double_t>   y(*data,"row_wise_branch.y");
    TTreeReaderValue<Int_t>   wom(*data,"row_wise_branch.WOMnumber");
    TTreeReaderValue<Double_t>   time(*data,"row_wise_branch.time");
    TTreeReaderValue<Int_t>   eventN(*data,"row_wise_branch.eventNumber");

    int entries = chain->GetEntries();
    int events = 10;

    TH2D *h1 = new TH2D("", "", 1000, -400, 400, 1000, -600, 600);
    h1->GetXaxis()->SetTitle("X (mm)");
    h1->GetYaxis()->SetTitle("Y (mm)");

    TH1I *h2 = new TH1I("", "", 2, 1, 3);
    h2->GetXaxis()->SetTitle("WOM Number");
    h2->GetYaxis()->SetTitle("Events");

    TH1D *h3 = new TH1D("", "", 100, 0, 1000);
    h3->GetXaxis()->SetTitle("Time (ns)");
    h3->GetYaxis()->SetTitle("Events");

    double N[events], w1[events], w2[events];
    for(int i = 0; i < events; i++) N[i] = i+1;

    while (data->Next()) {
      h1->Fill(*x,*y);
      h2->Fill(*wom);
      h3->Fill(*time);
      if(*wom == 1) w1[*eventN]++;
      if(*wom == 2) w2[*eventN]++;
    }

    ifstream fin("b/run1.mac");
    string s;
    s.reserve(50);
    string pos[2], dir[2];
    for(int i = 0; i < 2; ++i) std::getline(fin, s);
    std::getline(fin,s);
    int j = 14;
    for(int i = j; i < s.size(); ++i) {
      j++;
      if(s[i] == ' ') break;
      pos[0].push_back(s[i]);
    }
    for(int i = j; i < s.size(); ++i) {
      j++;
      if(s[i] == ' ') break;
      pos[1].push_back(s[i]);
    }
    std::getline(fin,s);
    j = 15;
    for(int i = j; i < s.size(); ++i) {
      j++;
      if(s[i] == ' ') break;
      dir[0].push_back(s[i]);
    }
    for(int i = j; i < s.size(); ++i) {
      j++;
      if(s[i] == ' ') break;
      dir[1].push_back(s[i]);
    }

    TCanvas *c1 = new TCanvas("c1","c1",10,10,800,1000);
    TLine* track = new TLine(stoi(pos[0]), stoi(pos[1]), stoi(pos[0]) + stoi(dir[0]) * 5000, stoi(pos[1]) + stoi(dir[1]) * 5000);
    track->SetLineColor(kRed);
    c1->SetLeftMargin(0.15);
    h1->SetStats(0);
    h1->Draw("COLZ");
    track->Draw("SAME");

    TCanvas *c2 = new TCanvas("c2","c2",10,10,1200,800);
    //gPad->SetLogy();
    h2->SetLineColor(kBlue);
    h2->Draw("HIST");

    TCanvas *c3 = new TCanvas("c3","c3",10,10,1200,800);
    gPad->SetLogy();
    h3->SetLineColor(kBlue);
    h3->Draw("HIST");

    TCanvas *c4 = new TCanvas("c4","c4",10,10,1200,800);
    c4->Divide(1,2);
    c4->cd(1);
    TGraph *h4 = new TGraph(events, N, w2);
    h4->GetXaxis()->SetTitle("Event");
    h4->GetYaxis()->SetTitle("Detected photons in Upper WOM");
    //gPad->SetLogy();
    h4->SetFillColor(38);
    h4->Draw("AB");
    c4->cd(2);
    TGraph *h5 = new TGraph(events, N, w1);
    h5->GetXaxis()->SetTitle("Event");
    h5->GetYaxis()->SetTitle("Detected photons in Lower WOM");
    //gPad->SetLogy();
    h5->SetFillColor(38);
    h5->Draw("AB");

    return;
}

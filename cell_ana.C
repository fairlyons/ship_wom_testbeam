void cell_ana() {
    gStyle->SetOptStat("eou");

    TString f = "b/cerenkov.root";
    TChain *chain = new TChain("Detected");
    chain->Add(f);

    TTreeReader *data = new TTreeReader(chain);
    TTreeReaderValue<Double_t>   x(*data,"row_wise_branch.x");
    TTreeReaderValue<Double_t>   y(*data,"row_wise_branch.y");
    TTreeReaderValue<Int_t>   WOM(*data,"row_wise_branch.WOMnumber");
    TTreeReaderValue<Double_t>   time(*data,"row_wise_branch.time");

    TH2D *h1 = new TH2D("", "", 1000, -400, 400, 1000, -600, 600);
    h1->GetXaxis()->SetTitle("X (mm)");
    h1->GetYaxis()->SetTitle("Y (mm)");

    TH1I *h2 = new TH1I("", "", 2, 1, 3);
    h2->GetXaxis()->SetTitle("WOM Number");
    h2->GetYaxis()->SetTitle("Events");

    TH1D *h3 = new TH1D("", "", 100, 0, 1000);
    h3->GetXaxis()->SetTitle("Time (?)");
    h3->GetYaxis()->SetTitle("Events");

    while (data->Next()) {
	h1->Fill(*x,*y);
        h2->Fill(*WOM);
        h3->Fill(*time);
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

    return;
}

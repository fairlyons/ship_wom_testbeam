#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

int cell_ana(){

    TCanvas* c1 = new TCanvas();
    TFile* f = new TFile("../build_ship_wom/test.root");
    TTree* t = (TTree*)f->Get("Detected");

    t->Draw("y:x>>h(600,-1200,1200,400,-800,800)","","colz");

    TLine* vline = new TLine(0,-800,0,800);
    TLine* hline = new TLine(-1200,0,1200,0);
    TLine* track = new TLine(300, 900, 300 - 1600*0.65, -800);
    track->SetLineColor(kRed);
    hline->Draw(); vline->Draw(); track->Draw();

    TCanvas* c2 = new TCanvas();
    c2->Divide(2,2);
    c2->cd(1);
    t->Draw("time>>h0(100,0,60)","WOWnumber==1");
    c2->cd(2);
    t->Draw("time>>h1(100,0,60)","WOWnumber==0");
    c2->cd(3);
    t->Draw("time>>h2(100,0,60)","WOWnumber==2");
    c2->cd(4);
    t->Draw("time>>h3(100,0,60)","WOWnumber==3");

    // TH1F *h1 = (TH1F*)c2->GetPrimitive("h1");
    // std::cout << h1->GetEntries();
    TCanvas* c3 = new TCanvas();
    c3->cd();
    TH1I* hist = new TH1I("hist", "Number of photons in each cell", 4, 0, 4);
    for(int i = 0; i < 120; i++)
        hist->Fill(0);

    for(int i = 0; i < 608; i++)
        hist->Fill(1);

    for(int i = 0; i < 1684; i++)
        hist->Fill(2);

    for(int i = 0; i < 0; i++)
        hist->Fill(3);


    hist->Draw();
    return 0;
}
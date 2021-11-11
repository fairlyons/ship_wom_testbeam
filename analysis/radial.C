#include "TTree.h"
#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"

void radial(){
    gROOT->ProcessLine(".L OrlovStyle.C");
    // gStyle->SetFrameFillColor(1);

    TFile *f = new TFile("../../build-SHiP_SBT_LScin/angular_dist_large.root");
    TTree *t = (TTree*)f->Get("Detected");

    TH1D *radial = new TH1D("radial",  "Radial distribution of photon hits", 
                            200, 26, 31);

    Double_t x, y;

    Int_t nentries = (Int_t)t->GetEntries();

    for (Int_t i=0; i<nentries; i++) {
        t->GetEntry(i);
        x = t->GetLeaf("x")->GetValue(0);
        y = t->GetLeaf("y")->GetValue(0);
        (x > 0) ? radial->Fill(sqrt((x-400)*(x-400) + y*y)) : 
                  radial->Fill(sqrt((x+400)*(x+400) + y*y));
    }

    
    TCanvas *myc = new TCanvas("myc", "myc");
    myc->SetLogy();
    radial->Draw();
    // myc->SaveAs("tube_hits.tex");
}
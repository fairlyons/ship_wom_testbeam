#include "TFile.h"
#include "TLine.h"
#include "TTree.h"

int cell_ana(){

    TFile* f = new TFile("../build_ship_wom/testref4.root");
    TTree* t = (TTree*)f->Get("Detected");

    t->Draw("y:x>>h(1200,-1200,1200,800,-800,800)","","colz");

    TLine* vline = new TLine(0,-800,0,800);
    TLine* hline = new TLine(-1200,0,1200,0);
    TLine* track = new TLine(1200, 510, 1300 + (1./0.7)*(-800-510), -800);
    track->SetLineColor(kRed);
    hline->Draw(); vline->Draw(); track->Draw();

    return 0;
}
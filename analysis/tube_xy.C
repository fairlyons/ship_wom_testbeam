#include "TTree.h"
#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"

void tube_xy(){
    gROOT->ProcessLine(".L OrlovStyle.C");

  // histogram divisions: only 5 in x to avoid label overlaps
    gStyle->SetNdivisions(10,"x");
    gStyle->SetNdivisions(10,"y");
    // gStyle->SetPalette(53);

    TCanvas *myc = new TCanvas("myc", "myc");


    // myc->Divide(3, 1, 0.001, 0.001);
    myc->SetWindowSize(1500,1000);

    TFile *f = new TFile("../../build-SHiP_SBT_LScin/optical_border.root");
    TTree *t = (TTree*)f->Get("Detected");
 
    Int_t nBins = 200;
    TH2D *ltube = new TH2D("left",  "Hits distribution at the left tube", 
                            nBins, -40, 40, nBins, -40,40);

    TH2D *rtube = new TH2D("right",  "Hits distribution at the right tube", 
                            nBins, -40, 40, nBins, -40,40);

    TH1D *radial = new TH1D("radial",  "Radial distribution of photon hits", 
                            nBins, 26, 31);

    TH1D *phi = new TH1D("phi",  "Phi distribution of photon hits", 
                            100, -180, 180);

    Double_t x, y;

    Int_t nentries = (Int_t)t->GetEntries();

    Double_t dist = 300;

    for (Int_t i=0; i<nentries; i++) {
        t->GetEntry(i);
        x = t->GetLeaf("x")->GetValue(0);
        y = t->GetLeaf("y")->GetValue(0);
        radial->Fill(x*x+y*y);
        (x > 0) ? rtube->Fill(x-dist,y) : ltube->Fill(x+dist,y);
        (x > 0) ? radial->Fill(sqrt((x-dist)*(x-dist) + y*y)) : 
                  radial->Fill(sqrt((x+dist)*(x+dist) + y*y));
        if (x > 0){
            x -= dist;
            if (x > 0)
                phi->Fill(atan(y/x) * 180. / 3.1415 );
            if (x < 0 && y > 0)
                phi->Fill(180. - atan(-y/x) * 180. / 3.1415 );
            if (x < 0 && y < 0)
                phi->Fill(-180. + atan(y/x) * 180. / 3.1415 );
        } else {
            x += dist;
            if (x > 0)
                phi->Fill(atan(y/x) * 180. / 3.1415 );
            if (x < 0 && y > 0)
                phi->Fill(180. - atan(-y/x) * 180. / 3.1415 );
            if (x < 0 && y < 0)
                phi->Fill(-180. + atan(y/x) * 180. / 3.1415 );
        }

    }



    TPad* leftPad = new TPad("leftPad", "leftPad", 
          .005, .505, .495, .995);
    leftPad->Draw();

    TPad* rightPad = new TPad("rightPad", "rightPad", 
          .505, .505, .995, .995);
    rightPad->Draw();

    TPad* lowPad = new TPad("lowPad", "lowPad", 
          .005, .005, .495, .45);
    lowPad->Draw();

    TPad* lowRightPad = new TPad("lowRightPad", "lowRightPad", 
          .505, .005, .995, .45);
    lowRightPad->Draw();

    leftPad->cd();
    leftPad->SetFrameFillColor(1);
    Double_t lnorm = 1/ltube->GetMaximum();
    ltube->Scale(lnorm);
    ltube->SetAxisColor(kGray+3, "XY");
    ltube->GetXaxis()->SetTitle("x, mm");
    ltube->GetYaxis()->SetTitle("y, mm");
    ltube->Draw("colz");

    rightPad->cd();
    rightPad->SetFrameFillColor(1);
    Double_t rnorm = 1/rtube->GetMaximum();
    rtube->Scale(rnorm);
    rtube->SetAxisColor(kGray+3, "XY");
    rtube->GetXaxis()->SetTitle("x, mm");
    rtube->GetYaxis()->SetTitle("y, mm");
    rtube->Draw("colz");


    lowPad->cd();
	radial->SetNdivisions(510,"X");
    Double_t lownorm = radial->Integral();
    // radial->Scale(1./lownorm); 
    radial->GetXaxis()->SetTitle("Radial distance, mm");
    radial->GetYaxis()->SetTitle("Counts");
    lowPad->SetLogy();
    radial->Draw();
	myc->Update();

	TLine *l1=new TLine(27.1,pow(10, lowPad->GetUymin()),27.1,pow(10, lowPad->GetUymax()));
	l1->SetLineColor(kBlack);
	l1->Draw();

	TLine *l2=new TLine(29.9,pow(10, lowPad->GetUymin()),29.9,pow(10, lowPad->GetUymax()));
	l2->SetLineColor(kBlack);
	l2->Draw();

	Double_t left_integ   = 100*radial->Integral(radial->GetXaxis()->FindBin(26), radial->GetXaxis()->FindBin(27.1))/lownorm   ;
	Double_t middle_integ = 100*radial->Integral(radial->GetXaxis()->FindBin(27.1), radial->GetXaxis()->FindBin(29.9))/lownorm ; 
	Double_t right_integ  = 100*radial->Integral(radial->GetXaxis()->FindBin(29.9), radial->GetXaxis()->FindBin(30.5))/lownorm ;

	TLatex l;
    l.DrawLatex(26.25, 1.1e4, Form("%.2f%%", left_integ));
	l.DrawLatex(28.25, 1.1e4, Form("%.2f%%", middle_integ));
    l.DrawLatex(30.1, 1.1e4, Form("%.2f%%", right_integ));


    lowRightPad->cd();
	phi->SetNdivisions(510,"X");
    phi->GetXaxis()->SetTitle("Phi, deg");
    phi->GetYaxis()->SetTitle("Counts");
    // lowRightPad->SetLogy();
    phi->Draw();
	myc->Update();



    myc->SaveAs("optical_border.tex");
}
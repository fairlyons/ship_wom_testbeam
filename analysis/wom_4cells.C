#include "TTree.h"
#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include <vector>
#include <iostream>
void wom_4cells()
{
  gROOT->ProcessLine(".L OrlovStyle.C");

  // histogram divisions: only 5 in x to avoid label overlaps
  gStyle->SetNdivisions(10, "x");
  gStyle->SetNdivisions(10, "y");
  // gStyle->SetPalette(53);

  TCanvas *myc = new TCanvas("myc", "myc");

  // myc->Divide(3, 1, 0.001, 0.001);
  // myc->SetWindowSize(1800, 1000);

  TFile *f = new TFile("../../build-SHiP_SBT_LScin/cerenkov.root");
  TTree *t = (TTree *)f->Get("Detected");

  std::vector<TPad *> pads;
  std::vector<TH2D *> histos;

  char name[20];

  int index = 0;
  for (int j = 0; j < 2; j++)
  {
    for (int i = 0; i < 4; i++)
    {
      sprintf(name, "%s%d","histo_",index);
      int nbins = 150;
      histos.push_back( new TH2D(name, name, nbins, -900 + i * 600 - 30, -900 + i * 600 + 30, 
                                             nbins, 400 - 800 * j - 30, 400 - 800 * j + 30) );
      index++;
    }
  }

  index = 0;

  myc->Divide(4, 2, 0.02, 0.02);
  for (int j = 0; j < 2; j++)
  {
    for (int i = 0; i < 4; i++)
    {
      sprintf(name, "%s%d","y:x>>histo_",index);
      myc->cd(index + 1);
      t->Draw(name, "", "colz");
      // histos[index]->Draw();
      myc->Update();

      index++;
    }
  }

  // myc->cd(1);
  // t->Draw(name, "", "colz");

  myc->Draw();
  // myc->SaveAs("optical_border.tex");
}
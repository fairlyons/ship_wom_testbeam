#include "TTree.h"
#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TCut.h"
#include <vector>
#include <iostream>
void wom_4cells(TString file_name)
{
  gROOT->ProcessLine(".L OrlovStyle.C");

  // histogram divisions: only 5 in x to avoid label overlaps
  gStyle->SetNdivisions(10, "x");
  gStyle->SetNdivisions(10, "y");
  // gStyle->SetPalette(53);

  TCanvas *myc = new TCanvas("myc", "myc");

  // myc->Divide(3, 1, 0.001, 0.001);
  myc->SetWindowSize(1500, 750);

  TFile *f = new TFile(file_name);
  TTree *t = (TTree *)f->Get("Detected");

  std::vector<TPad *> pads;
  std::vector<TH2D *> histos;
  std::vector<TH1D *> time_histos;

  char name[20];
  char title[20];
  int index = 0;
  for (int j = 0; j < 2; j++)
  {
    for (int i = 0; i < 4; i++)
    {
      sprintf(name, "%s%d","histo_",index);
      sprintf(title, "%s%d","WOM ",index);
      int nbins = 150;
      histos.push_back( new TH2D(name, title, nbins, -900 + i * 600 - 30, -900 + i * 600 + 30, 
                                             nbins, 400 - 800 * j - 30, 400 - 800 * j + 30) );
     
      sprintf(name, "%s%d","time_histo_", index);
      sprintf(title, "%s%d","Time in WOM ", index);

      time_histos.push_back( new TH1D(name, title, 300, 0, 40) );
      index++;
    }
  }

  index = 0;

  myc->Divide(4, 2, 0.02, 0.02);
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(4, 2, 0.02, 0.02);


  TCut rightWOM = "sipmNumber > 39";
  TCut leftWOM = "sipmNumber < 40";
  TCut detection = "detection == 1";

  TCut cell_0 = "WOMnumber == 0";
  TCut cell_1 = "WOMnumber == 1";
  TCut cell_2 = "WOMnumber == 2";
  TCut cell_3 = "WOMnumber == 3";

  std::vector<TCut> wom_cut;

  wom_cut.push_back(cell_1 &&  leftWOM);
  wom_cut.push_back(cell_1 && rightWOM);
  wom_cut.push_back(cell_0 &&  leftWOM);
  wom_cut.push_back(cell_0 && rightWOM);
  wom_cut.push_back(cell_2 &&  leftWOM);
  wom_cut.push_back(cell_2 && rightWOM);
  wom_cut.push_back(cell_3 &&  leftWOM);
  wom_cut.push_back(cell_3 && rightWOM);

  for (int j = 0; j < 2; j++)
  {
    for (int i = 0; i < 4; i++)
    {
      sprintf(name, "%s%d","y:x>>histo_",index);
      myc->cd(index + 1);
      t->Draw(name, "detection==1", "colz");
      myc->Update();

      sprintf(name, "%s%d","time>>time_histo_",index);
      c2->cd(index + 1);
      TCut cut = wom_cut[index] && detection;
      t->Draw(name, cut,  "same");
      c2->Update();

      index++;
    }
  }

  // myc->cd(1);
  // t->Draw(name, "", "colz");

  myc->Draw();




  // myc->SaveAs("optical_border.tex");
}
void An_photons() {
    
    TFile *outfile  = TFile::Open("PH.root","RECREATE");
    int Category;
    TTree *tree = new TTree("T","CERN 1988 staff data");
    tree->Branch("x-y",&Category,"Category/I");
    tree->Branch("WOM",&Category,"Category/I");
    tree->Branch("time",&Category,"Category/I");
for(int G=0;  G<2;  G++){
    string GG = to_string(G);
    string name = "/home/marte/software/ship_wom_testbeam/bb/photons";
    name+=GG.c_str();
    name+=".root";
    TFile *f = new TFile(name.c_str());
    TTree *t = (TTree*)f->Get("Detected");
    
    double x, y, WOM, time, eventNumber;

    int nentries =t->GetEntries();

    TH2D *h1 = new TH2D("", "", 100, -400, 400, 100, -600, 600);
    h1->GetXaxis()->SetTitle("X (mm)");
    h1->GetYaxis()->SetTitle("Y (mm)");
    
    TH1I *h2 = new TH1I("", "", 2, 1, 3);
    h2->GetXaxis()->SetTitle("WOM Number");
    h2->GetYaxis()->SetTitle("Events");

    TH1D *h3 = new TH1D("", "", 100, 0, 100);
    h3->GetXaxis()->SetTitle("Time (ns)");
    h3->GetYaxis()->SetTitle("Events");
 //   double N[nentries], w1[nentries], w2[nentries];
  //  for(int i = 0; i < nentries; i++) N[i] = i+1;


    for (Int_t i=0; i<nentries; i++) {
        t->GetEntry(i);
        x = t->GetLeaf("x")->GetValue(0);
        y = t->GetLeaf("y")->GetValue(0);
        WOM = t->GetLeaf("WOMnumber")->GetValue(0);
        time = t->GetLeaf("time")->GetValue(0);
        eventNumber = t->GetLeaf("eventNumber")->GetValue(0);
        h1->Fill(x,y);
        h2->Fill(WOM);
        h3->Fill(time);
      //  if(WOM == 1) w1[i]++;
      //  if(WOM == 2) w2[i]++;
    }
    string name_c1 = "x-y";
    name_c1 += GG.c_str(); 
    TCanvas *c1 = new TCanvas(name_c1.c_str(),"c1",10,10,800,1000);
    c1->SetLeftMargin(0.15);
    h1->SetStats(0);
    h1->Draw("COLZ");
    string image = "Image/";
    image += name_c1; 
    image += ".png";
    c1->SaveAs(image.c_str());
    outfile->WriteObject(h1, "x-y");	
  
    string name_c2 = "WOMNumber";
    name_c2 += GG.c_str(); 
    TCanvas *c2 = new TCanvas(name_c2.c_str(),"c2",10,10,1200,800);
    h2->SetLineColor(kBlue);
    h2->Draw("HIST");
    image = "Image/";
    image += name_c2; 
    image += ".png";
    c2->SaveAs(image.c_str());		
    outfile->WriteObject(h2, "WOM");	

    string name_c3 = "time";
    name_c3 += GG.c_str(); 
    TCanvas *c3 = new TCanvas(name_c3.c_str(),"c3",10,10,1200,800);
    gPad->SetLogy();
    h3->SetLineColor(kBlue);
    h3->Draw("HIST");
    image = "Image/";
    image += name_c3; 
    image += ".png";
    c3->SaveAs(image.c_str());		
    outfile->WriteObject(h2, "time");		

    /* TCanvas *c4 = new TCanvas("c4","c4",10,10,1200,800);
    c4->Divide(1,2);
    c4->cd(1);
    TGraph *h4 = new TGraph(eventNumber, N, w2);
    h4->GetXaxis()->SetTitle("Event");
    h4->GetYaxis()->SetTitle("Detected photons in Upper WOM");
    //gPad->SetLogy();
    h4->SetFillColor(38);
    h4->Draw("AB");
    c4->cd(2);
   TGraph *h5 = new TGraph(eventNumber, N, w1);
    h5->GetXaxis()->SetTitle("Event");
    h5->GetYaxis()->SetTitle("Detected photons in Lower WOM");
    //gPad->SetLogy();
    h5->SetFillColor(38);
    h5->Draw("AB");*/
}
    return;
}



void Q2linearity(){
  
  //read in the text file with Iset and B[kG]
  TTree *tt = new TTree("tt", "tt");
  tt->ReadFile("Q2input.txt", "II:BB");
  int nentries = tt->GetEntries();
  cout<<"number of entries:\t"<<nentries<<endl;
  const int num = 14;
  double Iset[num];
  double Bfield[num];
  double BfieldE[num];
  
  for (int ii=0; ii<num; ii++){
      tt->GetEntry(ii);
      TLeaf *BfieldL = tt->GetLeaf("BB");
      TLeaf *IsetL = tt->GetLeaf("II");
      Bfield[ii] = BfieldL->GetValue();
      Iset[ii] = IsetL->GetValue();
      BfieldE[ii] = 0.0001;
      cout<<"entry:\t"<<ii<<"\t"<<Iset[ii]<<"\t"<<Bfield[ii]<<endl;
    }
  
  
  //make the output root file
  TFile *fout = new TFile("Q2LinearityFit.root","RECREATE");
  
  //plot B vs I with errors
  TGraphErrors *g1 = new TGraphErrors(num,Iset,Bfield,0,BfieldE);
  TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  g1->SetMarkerStyle(21);
  g1->SetMarkerColor(kBlue);
  g1->Fit(f1);
  g1->Draw();
  g1->Write("DataFit");

  //Use the fit to plot the absolute residuals
  double AbsRes[num];
  double RelRes[num];
  for (int ii=0; ii<num; ii++){
    AbsRes[ii] = Bfield[ii] - f1->Eval(Iset[ii]);
    RelRes[ii] = AbsRes[ii]/Bfield[ii];
  }
  TGraphErrors *g2 = new TGraphErrors(num,Iset,AbsRes,0,BfieldE);
  g2->SetMarkerStyle(22);
  g2->SetMarkerColor(kRed);
  g2->Draw();
  g2->Write("AbsRes");

  
  //Plot the relative residuals
  TGraphErrors *g3 = new TGraphErrors(num,Iset,RelRes,0,BfieldE);
  g3->SetMarkerStyle(23);
  g3->SetMarkerColor(kMagenta);
  g3->Draw();
  g3->Write("RelRes");

  fout->Write();
  fout->Close();




}

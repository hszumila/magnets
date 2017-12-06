double correctNegPol(double mkG){

  double corr = 0;
  return mkG + corr;

}

double correctPosPol(double mkG){

  double corr = 0;
  return mkG - corr;

}

double fasdf(double *xx,double *pp)
{
  //none work:
  //return (xx[0]<pp[0])*pp[1] + (xx[0]>pp[0])*(pp[2]*xx[0]*xx[0]+pp[3]*xx[0]+pp[4]);
  // return (xx[0]<=pp[0])*pp[1] + (xx[0]>pp[0])*(pp[1]-(pp[2]*xx[0]*xx[0]+pp[3]*xx[0]+pp[4]));
  //return pp[0] - pp[1]*exp(-(xx[0]-pp[2])*pp[3]);
  //return pp[0]*TMath::Erfc(pp[1]*(xx[0]+pp[2]));

  //for 2nd order poly
  //return pp[1] + (xx[0]>pp[0])*pp[2]*pow(xx[0]-pp[0],2);

  //for gaussian tail
  //return pp[1]*(xx[0]<pp[0]) + (xx[0]>pp[0])*exp(-0.5*pow((xx[0]-pp[0])/pp[2],2));

   //for high poly tail, set par = 5
  return pp[1] + (xx[0]>pp[0])*(pp[2]*pow(xx[0]-pp[0],3)+pp[3]*pow(xx[0]-pp[0],2)+pp[4]*(xx[0]-pp[0]));
}

double fasdf2(double *xx,double *pp)
{
 
  return (xx[0]<=pp[0])*(pp[1]+pp[2]*xx[0]+pp[3]*pow(xx[0],2)+pp[4]*pow(xx[0],3))+(xx[0]>pp[0])*(pp[5]+pp[6]*(xx[0]-pp[0])+pp[7]*pow(xx[0]-pp[0],2)+pp[8]*pow(xx[0]-pp[0],3));
   //for high poly tail, set par = 5
  // return pp[1] + (xx[0]>pp[0])*(pp[2]*pow(xx[0]-pp[0],3)+pp[3]*pow(xx[0]-pp[0],2)+pp[4]*(xx[0]-pp[0])) + (xx[0]<pp[0])*(pp[5]*xx[0]);
}

void plotDataExt(){

  //read the input: +ramp1, +ramp2, -ramp, error table
  TTree *r1 = new TTree("r1","r1");
  r1->ReadFile("input/ramp1ext.txt","c1:ramp1");
  const int nRamp1Pts = r1->GetEntries();
  double r1_B[nRamp1Pts];
  double r1_I[nRamp1Pts];
  for (int jj=0;jj<nRamp1Pts; jj++){
    r1->GetEntry(jj);
    TLeaf *r1BL = r1->GetLeaf("ramp1");
    TLeaf *r1IL = r1->GetLeaf("c1");
    r1_B[jj] = r1BL->GetValue();//kG
    r1_I[jj] = r1IL->GetValue();
  }

  TTree *r2 = new TTree("r2","r2");
  r2->ReadFile("input/ramp2ext.txt","c2:ramp2");
  const int nRamp2Pts = r2->GetEntries();
  double r2_B[nRamp2Pts];
  double r2_I[nRamp2Pts];
  for (int jj=0;jj<nRamp2Pts; jj++){
    r2->GetEntry(jj);
    TLeaf *r2BL = r2->GetLeaf("ramp2");
    TLeaf *r2IL = r2->GetLeaf("c2");
    r2_B[jj] = r2BL->GetValue();//kG
    r2_I[jj] = r2IL->GetValue();
  }

  TTree *r = new TTree("r","r");
  r->ReadFile("input/rampext.txt","c:ramp");
  const int nRampPts = r->GetEntries();
  double r_B[nRampPts];
  double r_I[nRampPts];
  double r_Ip[nRampPts];
  for (int jj=0;jj<nRampPts; jj++){
    r->GetEntry(jj);
    TLeaf *rBL = r->GetLeaf("ramp");
    TLeaf *rIL = r->GetLeaf("c");
    r_B[jj] = -rBL->GetValue();//kG
    r_I[jj] = -rIL->GetValue();
    r_Ip[jj] = rIL->GetValue();
  }

  TTree *r3 = new TTree("r3","r3");
  r3->ReadFile("input/rampCombext.txt","c3:ramp3");
  const int nRampPts3 = r3->GetEntries();
  double r3_B[nRampPts3];
  double r3_I[nRampPts3];
  for (int jj=0;jj<nRampPts3; jj++){
    r3->GetEntry(jj);
    TLeaf *rBL3 = r3->GetLeaf("ramp3");
    TLeaf *rIL3 = r3->GetLeaf("c3");
    r3_B[jj] = rBL3->GetValue();//kG
    r3_I[jj] = rIL3->GetValue();
    cout<<"jj:\t"<<jj<<"\tcurrent:\t"<<r3_I[jj]<<"\tfield:\t"<<r3_B[jj]<<endl;

    
  }

  
  TTree *l = new TTree("l","l");
  l->ReadFile("input/leff.txt","cl:leff");
  const int ncl = l->GetEntries();
  double l_leff[ncl];
  double l_c[ncl];
  for (int jj=0;jj<ncl; jj++){
    l->GetEntry(jj);
    TLeaf *leffL = l->GetLeaf("leff");
    TLeaf *lcL = l->GetLeaf("cl");
    l_c[jj] = lcL->GetValue();
    l_leff[jj] = leffL->GetValue();
  }

   TTree *bl = new TTree("bl","bl");
  bl->ReadFile("input/oldBetaLeff.txt","cc:bl");
  const int nbl = bl->GetEntries();
  double l_bleff[nbl];
  double l_cc[nbl];
  for (int jj=0;jj<nbl; jj++){
    bl->GetEntry(jj);
    TLeaf *bleffL = bl->GetLeaf("bl");
    TLeaf *ccL = bl->GetLeaf("cc");
    l_cc[jj] = ccL->GetValue();
    l_bleff[jj] = bleffL->GetValue();
  }

  
  //make output pdf
  TFile *fout = new TFile("output/hms_dipExt.root","RECREATE");
  TCanvas *canvas = new TCanvas("canvas","hms dipole data", 800, 800);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin(0.2);
  //canvas->SetRightMargin(0.2);
  std::string pdf_file_name = "output/hms_dipExt.pdf";
  gROOT->SetBatch(true);

  canvas->Update();

  
  //calc the hall probe offset and remant field
  double probeOffset = 0.0;//(r1_B[nRamp1Pts-1]+r_B[nRampPts-1])/2;
  double remField = 0.0;//(r1_B[nRamp1Pts-1]-r_B[nRampPts-1])/2;

  cout<<"Probe offset:\t"<<probeOffset<<" kG\t remnant field:\t"<<remField<<" kG"<<endl;
  
  //correct the data for the probe offset and error correction
  double corr_Bpos[nRamp1Pts];
  double corr_Bneg[nRampPts];
  double asyCorr[nRampPts];
  double asyNoCorr[nRampPts];
  double asyCorrRel[nRampPts];
  double asyNoCorrRel[nRampPts];
  for (int ii=0; ii<nRamp1Pts;ii++){
    corr_Bpos[ii] = correctPosPol(r1_B[ii]);
    if (r1_I[ii]<1000){corr_Bpos[ii] = correctPosPol(r1_B[ii] - probeOffset);}
    
    for (int jj=0; jj<nRampPts; jj++){
      if (abs(r_I[jj]) == r1_I[ii]){
	corr_Bneg[jj] = correctNegPol(r_B[jj]);
	asyNoCorr[jj] = r_B[jj] + r1_B[ii];

	if (abs(r_I[jj])<1000){
	  corr_Bneg[jj] = correctNegPol(r_B[jj] - probeOffset);
	  asyNoCorr[jj] = r_B[jj]-probeOffset + r1_B[ii] - probeOffset;
	}
	
	asyCorr[jj] = corr_Bneg[jj] + corr_Bpos[ii];
	asyNoCorrRel[jj] = asyNoCorr[jj]/corr_Bpos[ii];
	asyCorrRel[jj] = asyCorr[jj]/corr_Bpos[ii];
      }
    }
  }  

  double BmaxPos = corr_Bpos[0];
  double BmaxNeg = corr_Bneg[0];

  //plot the residual with straight line
  double resPos[nRamp1Pts];
  double relResPos[nRamp1Pts];
  for (int ii=0; ii<nRamp1Pts; ii++){
    resPos[ii] = corr_Bpos[ii] - BmaxPos*r1_I[ii]/r1_I[0];
    relResPos[ii] = resPos[ii]/corr_Bpos[ii];
  }
  double resNeg[nRampPts];
  double relResNeg[nRampPts];
  for (int ii=0; ii<nRampPts;ii++){
    resNeg[ii] = corr_Bneg[ii] - BmaxNeg*r_I[ii]/r_I[0];
    relResNeg[ii] = resNeg[ii]/corr_Bneg[ii];
  }

  //plot the hysteresis
  TGraph *g1 = new TGraph(nRamp1Pts,r1_I,r1_B);
  g1->SetMarkerStyle(22);
  g1->SetMarkerColor(kBlue);
  g1->SetMarkerSize(1.5);
  TText *t_g1 = new TText(500,-5,"pos ramp 1");
  t_g1->SetTextColor(kBlue);

  TGraph *g2 = new TGraph(nRamp2Pts,r2_I,r2_B);
  g2->SetMarkerStyle(22);
  g2->SetMarkerColor(kRed);
  g2->SetMarkerSize(1.5);
  TText *t_g2 = new TText(500,-10,"pos ramp 2");
  t_g2->SetTextColor(kRed);

  TGraph *g3 = new TGraph(nRampPts,r_I,r_B);
  g3->SetMarkerStyle(22);
  g3->SetMarkerColor(kCyan+1);
  g3->SetMarkerSize(1.5);
  g3->SetTitle("B vs I");
  g3->GetXaxis()->SetTitle("Current [A]");
  g3->GetYaxis()->SetTitle("measured B field [kG]");
  g3->GetYaxis()->SetTitleOffset(2.0);
  TText *t_g3 = new TText(500,-15,"neg ramp");
  t_g3->SetTextColor(kCyan+1);

  g3->Draw("ap");
  g3->GetXaxis()->SetLimits(-3100,3100);
  g3->SetMinimum(-24);
  g3->SetMaximum(24);
  g3->Draw("p");
  canvas->Update();
  g1->Draw("psame");
  g2->Draw("psame");
  TF1 *f1 = new TF1("f1","[0]+[1]*x",-2000,2000);
  f1->SetLineStyle(2);
  f1->SetLineColor(kMagenta);
  g1->Fit(f1,"QR");
  TText *t_line = new TText(-2100,15,Form("B [kG] = %.3f I + %.3f",f1->GetParameter(1),f1->GetParameter(0)));
  t_line->SetTextColor(kMagenta);
  t_line->Draw();
  t_g1->Draw();
  t_g2->Draw();
  t_g3->Draw();
  
  canvas->Print( (pdf_file_name + "(").c_str());

  //plot the rel residual to straight line
  TGraph *g4 = new TGraph(nRamp1Pts,r1_I,resPos);
  g4->SetMarkerStyle(22);
  g4->SetMarkerColor(kBlue);
  g4->SetMarkerSize(1.5);
  TGraph *g5 = new TGraph(nRampPts,r_I,resNeg);
  g5->SetMarkerStyle(22);
  g5->SetMarkerColor(kRed);
  g5->SetMarkerSize(1.5);

  g4->SetTitle("Residual relative to line from (0,0) to (B_{max}, I_{max})");
  g4->GetXaxis()->SetTitle("Current [A]");
  g4->GetYaxis()->SetTitle("Residual [kG]");
  g4->GetYaxis()->SetTitleOffset(1.8);
  g4->Draw("ap");
  g4->GetXaxis()->SetLimits(-3100,3100);
  g4->SetMinimum(-5);
  g4->SetMaximum(5);
  g5->Draw("psame");

  canvas->Print( (pdf_file_name + "(").c_str());

  TGraph *g6 = new TGraph(nRamp1Pts,r1_I,relResPos);
  g6->SetMarkerStyle(22);
  g6->SetMarkerColor(kBlue);
  g6->SetMarkerSize(1.5);
  TGraph *g7 = new TGraph(nRampPts,r_I,relResNeg);
  g7->SetMarkerStyle(22);
  g7->SetMarkerColor(kRed);
  g7->SetMarkerSize(1.5);

  g6->SetTitle("Relative Residual relative to line from (0,0) to (B_{max}, I_{max})");
  g6->GetXaxis()->SetTitle("Current [A]");
  g6->GetYaxis()->SetTitle("Relative Residual [unitless]");
  g6->GetYaxis()->SetTitleOffset(1.8);
  g6->Draw("ap");
  g6->GetXaxis()->SetLimits(-3100,3100);
  g6->SetMinimum(0);
  g6->SetMaximum(0.5);
  g7->Draw("psame");

  canvas->Print( (pdf_file_name + "(").c_str());
  
  //plot the asymmetry analysis with and without the probe error corrections
  TGraph *gasy = new TGraph(nRampPts,r_Ip,asyNoCorr);
  gasy->SetMarkerStyle(22);
  gasy->SetMarkerColor(kBlue);
  gasy->SetMarkerSize(1.5);
  TGraph *gasyc = new TGraph(nRampPts,r_Ip,asyCorr);
  gasyc->SetMarkerStyle(22);
  gasyc->SetMarkerColor(kRed);
  gasyc->SetMarkerSize(1.5);

  TText *t_asyn = new TText(100,0.1,"not probe error corrected");
  TText *t_asyc = new TText(100,0.08,"probe error corrected");
  t_asyn->SetTextColor(kBlue);
  t_asyc->SetTextColor(kRed);
  
  gasy->SetTitle("Polarity asymmetry");
  gasy->GetXaxis()->SetTitle("Current [A]");
  gasy->GetYaxis()->SetTitle("B_{pos} + B_{neg} [kG]");
  gasy->GetYaxis()->SetTitleOffset(1.8);
  gasy->Draw("ap");
  gasy->SetMinimum(-0.01);
  gasy->SetMaximum(0.01);
  gasyc->Draw("psame");
  // t_asyn->Draw();
  //t_asyc->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  //plot the relative polarity asymmetry
  TGraph *gasyrel = new TGraph(nRampPts,r_Ip,asyNoCorrRel);
  gasyrel->SetMarkerStyle(22);
  gasyrel->SetMarkerColor(kBlue);
  gasyrel->SetMarkerSize(1.5);
  TGraph *gasycrel = new TGraph(nRampPts,r_Ip,asyCorrRel);
  gasycrel->SetMarkerStyle(22);
  gasycrel->SetMarkerColor(kRed);
  gasycrel->SetMarkerSize(1.5);

  TText *t_asynr = new TText(500,0.003,"not probe error corrected");
  TText *t_asycr = new TText(500,0.0025,"probe error corrected");
  t_asynr->SetTextColor(kBlue);
  t_asycr->SetTextColor(kRed);
  
  gasyrel->SetTitle("Relative polarity asymmetry");
  gasyrel->GetXaxis()->SetTitle("Current [A]");
  gasyrel->GetYaxis()->SetTitle("B_{pos} + B_{neg} / B_{pos} [unitless]");
  gasyrel->GetYaxis()->SetTitleOffset(1.8);
  gasyrel->Draw("ap");
  gasyrel->SetMinimum(-0.001);
  gasyrel->SetMaximum(0.001);
  gasycrel->Draw("psame");
  //t_asynr->Draw();
  //t_asycr->Draw();
  canvas->Print( (pdf_file_name + "(").c_str());

  //plot B/I vs I and fit well, plot res
  //plot B/I vs I
  const double beta_g =  8.93508498183325471E-03;
  double error[nRamp1Pts];
  double BIratio[nRamp1Pts];
  double BIratio3[nRampPts3];
  double BIratioG[nRamp1Pts];
  for (int ii=0; ii<nRamp1Pts; ii++){
    BIratio[ii] = corr_Bpos[ii]/r1_I[ii];
    BIratioG[ii] = BIratio[ii]/beta_g;
    if (ii == nRamp1Pts-1){
      BIratio[ii] = BIratio[ii-1];
      BIratioG[ii] = BIratio[ii-1]/beta_g;
    }
    error[ii] =0.0;// 0.001*BIratio[ii];
  }

  for (int jj=0; jj<nRampPts3; jj++){
    cout<<"jj:\t"<<jj<<"\tcurrent:\t"<<r3_I[jj]<<"\tfield:\t"<<r3_B[jj]<<endl;
    BIratio3[jj] = r3_B[jj]/r3_I[jj];
    cout<<"\tratio:\t"<<BIratio3[jj]<<endl;


  }

  
  /*
  gStyle->SetOptFit(111);
  TGraphErrors *gBI = new TGraphErrors(nRamp1Pts,r1_I,BIratio,0,error);
  gBI->SetMarkerStyle(22);
  gBI->SetMarkerColor(kBlue);
  gBI->SetMarkerSize(1.5);
  gBI->SetTitle("B vs I");
  gBI->GetXaxis()->SetTitle("Current [A]");
  gBI->GetYaxis()->SetTitle("#beta, B/I [kG]");
  gBI->GetYaxis()->SetTitleOffset(2.2);
  gBI->Draw("ap");
  TF1 *fBI = new TF1("fBI","pol4",0,3000);
  gBI->Fit(fBI,"R");
  */

  //gStyle->SetOptFit(111);
  TGraphErrors *gBI = new TGraphErrors(nRampPts3,r3_I,BIratio3,0,0);
  gBI->SetMarkerStyle(22);
  gBI->SetMarkerColor(kBlue);
  gBI->SetMarkerSize(1.5);
  gBI->SetTitle("B vs I");
  gBI->GetXaxis()->SetTitle("Current [A]");
  gBI->GetYaxis()->SetTitle("#beta, B/I [kG]");
  gBI->GetYaxis()->SetTitleOffset(2.2);
  gBI->Draw("ap");
  TF1 *fBI = new TF1("fBI","pol8",100,2900);

  /*
 TF1 *fBI = new TF1("fBI",fasdf2,100,2900,9);
 fBI->SetParLimits(0,1200,1500);//1200-1300//1000-1300
 fBI->SetParameter(1,0.00908);
 fBI->SetParameter(2,-4.519E-7);
 fBI->SetParameter(3,5.204E-10);
 fBI->SetParameter(4,-2.138E-13);
 //fBI->SetParameter(5,0.0089);
 fBI->SetParLimits(5,0.00885,0.009);
 fBI->SetParameter(6,0.482E-6);
 fBI->SetParameter(5,-2.431E-9);
 fBI->SetParameter(8,3.358E-13);
  */
  /*
   TF1 *fBI = new TF1("fBI",fasdf2,100,2900,9);
 fBI->SetParLimits(0,1200,1500);//1200-1300//1000-1300
 fBI->SetParameter(1,0.00908);
 fBI->SetParameter(2,-4.519E-7);
 fBI->SetParameter(3,5.204E-10);
 fBI->SetParameter(4,-2.138E-13);
 //fBI->SetParameter(5,0.0089);
 fBI->SetParLimits(5,0.00886,0.009);
 fBI->SetParameter(6,1.49189E-7);
 fBI->SetParameter(5,-6.16454E-10);
 fBI->SetParameter(8,-1.46953E-13);
 
   
  
  gBI->Fit(fBI,"R0");
  TF1 *fBI2 = new TF1("fBI2",fasdf2,100,2900,9);
  fBI2->SetParameter(0,fBI->GetParameter(0));
 fBI2->SetParameter(1,fBI->GetParameter(1));
 fBI2->SetParameter(2,fBI->GetParameter(2));
 fBI2->SetParameter(3,fBI->GetParameter(3));
 fBI2->SetParameter(4,fBI->GetParameter(4));
 fBI2->SetParameter(5,fBI->GetParameter(5));
 fBI2->SetParameter(6,fBI->GetParameter(6));
 fBI2->SetParameter(7,fBI->GetParameter(7));
 fBI2->SetParameter(8,fBI->GetParameter(8));
  */
   gBI->Fit(fBI,"R");

  canvas->Print( (pdf_file_name + "(").c_str());
  gStyle->SetOptFit(0);

  double fitResidual[nRampPts3];
  double fitRelResidual[nRampPts3];
  for (int ii=0; ii<nRampPts3; ii++){
    fitResidual[ii] = BIratio3[ii] - fBI->Eval(r3_I[ii]);
    fitRelResidual[ii] = fitResidual[ii]/BIratio3[ii];
  }

  TGraph *gBIres = new TGraph(nRampPts3,r3_I,fitResidual);
  gBIres->SetMarkerStyle(22);
  gBIres->SetMarkerColor(kBlue);
  gBIres->SetMarkerSize(1.5);
  gBIres->SetTitle("B/I vs I, fit residual");
  gBIres->GetXaxis()->SetTitle("Current [A]");
  gBIres->GetYaxis()->SetTitle("#beta, B/I_{data} - B/I_{fit} [kG/A]");
  gBIres->GetYaxis()->SetTitleOffset(2.2);
  gBIres->Draw("ap");
  gBIres->SetMinimum(-0.00002);
  gBIres->SetMaximum(0.00002);
  canvas->Print( (pdf_file_name + "(").c_str());
  
 
  TGraph *gBIrelres = new TGraph(nRampPts3,r3_I,fitRelResidual);
  gBIrelres->SetMarkerStyle(22);
  gBIrelres->SetMarkerColor(kBlue);
  gBIrelres->SetMarkerSize(1.5);
  gBIrelres->SetTitle("B/I vs I, fit relative residual");
  gBIrelres->GetXaxis()->SetTitle("Current [A]");
  gBIrelres->GetYaxis()->SetTitle("B/I_{data} - B/I_{fit} / B/I_{data} [unitless]");
  gBIrelres->GetYaxis()->SetTitleOffset(2.2);
  gBIrelres->Draw("ap");
  gBIrelres->SetMinimum(-0.002);
  gBIrelres->SetMaximum(0.002);
  canvas->Print( (pdf_file_name + "(").c_str());
  
  
  //plot beta_measured/beta_golden vs I, define beta_golden
  TGraphErrors *betaRatioG = new TGraphErrors(nRamp1Pts, r1_I, BIratioG, 0, error);
  betaRatioG->SetMarkerStyle(22);
  betaRatioG->SetMarkerSize(1.5);
  betaRatioG->SetMarkerColor(kBlue);
  betaRatioG->SetTitle("B/I ratio to golden tune");
  betaRatioG->GetXaxis()->SetTitle("Current [A]");
  betaRatioG->GetYaxis()->SetTitle("#beta, B/I_{data} / B/I_{golden} [unitless]");
  betaRatioG->GetYaxis()->SetTitleOffset(2.4);
  betaRatioG->Draw("ap");
  TF1 *fBIrat = new TF1("fBIrat","pol6",100,3000);
  betaRatioG->Fit(fBIrat,"QR");
  

  canvas->Print( (pdf_file_name + "(").c_str());



  TGraph *gl = new TGraph (ncl,l_c,l_leff);
  gl->SetMarkerStyle(22);
  gl->SetMarkerSize(1.5);
  gl->SetMarkerColor(kBlue);
  gl->SetTitle("Ratio of leff_{TOSCA}/leff_{const}");
  gl->GetXaxis()->SetTitle("Current [A]");
  gl->GetYaxis()->SetTitle("leff_{TOSCA} / leff_{const} [unitless]");
  gl->GetYaxis()->SetTitleOffset(2.4);
  gl->Draw("ap");

  //2nd order poly:
  /*
  TF1 *fleff = new TF1("fleff",fasdf,0,2500,3);//4);
  fleff->SetParameter(0,1000);
  fleff->FixParameter(1,1);
  fleff->SetParameter(2,-5E-9);
  */
  
  //for gaussian tail:
  /*
  TF1 *fleff = new TF1("fleff",fasdf,0,2500,3);
  fleff->SetParLimits(0,900,1300);
  fleff->FixParameter(1,1);
  fleff->SetParameter(2,700);
  */

  TF1 *fleff = new TF1("fleff",fasdf,0,3000,5);
  fleff->FixParameter(0,1220);
  fleff->FixParameter(1,1);
  fleff->SetParameter(2,0.7965E-12);
  fleff->SetParameter(3,-1.21752E-8);
  fleff->SetParameter(4,1.90788E-6);
 
  
  gl->Fit(fleff,"R");
  canvas->Print( (pdf_file_name + "(").c_str());


  double lrelres[ncl];
  for (int ii=0; ii< ncl; ii++){
    lrelres[ii] = (l_leff[ii] - fleff->Eval(l_c[ii]))/l_leff[ii];
  }

  //plot leff, res
  TGraph *glr = new TGraph (ncl,l_c,lrelres);
  glr->SetMarkerStyle(22);
  glr->SetMarkerSize(1.5);
  glr->SetMarkerColor(kBlue);
  glr->SetTitle("Fit to leff, relative residual");
  glr->GetXaxis()->SetTitle("Current [A]");
  glr->GetYaxis()->SetTitle("leff_{ratio} - leff_{ratio,fit}  / leff_{ratio} [unitless]");
  glr->GetYaxis()->SetTitleOffset(2.4);
  glr->Draw("ap");
  canvas->Print( (pdf_file_name + "(").c_str());


  //multiple the beta times the leff ratio
  double BILeff[nRamp1Pts];
  for (int ii=0; ii<nRamp1Pts; ii++){
    BILeff[ii] = BIratio[ii]*fleff->Eval(r1_I[ii]);
  }


  //plot the old golden beta and the new beta times leff
  TGraph *gg1 = new TGraph(nRamp1Pts,r1_I,BILeff);
  gg1->SetMarkerStyle(22);
  gg1->SetMarkerColor(kBlue);
  gg1->SetMarkerSize(1.5);
  TText *t_gg1 = new TText(500,-5,"B/I * Leff, 2017");
  t_gg1->SetTextColor(kBlue);

  TGraph *gg2 = new TGraph(nbl,l_cc,l_bleff);
  gg2->SetMarkerStyle(22);
  gg2->SetMarkerColor(kRed);
  gg2->SetMarkerSize(1.5);
  TText *t_gg2 = new TText(500,-10,"B/I * Leff, OLD");
  t_gg2->SetTextColor(kRed);

  gg1->SetTitle("B/I * Leff ratio");
  gg1->GetXaxis()->SetTitle("Current [A]");
  gg1->GetYaxis()->SetTitle("B/I * L (ratio) [kG/A]");
  gg1->GetYaxis()->SetTitleOffset(1.8);
  gg1->Draw("ap");
  gg1->SetMinimum(0.006);
  gg1->SetMaximum(0.0095);
  gg2->Draw("psame");

  TF1 *fBILrat = new TF1("fBILrat","pol6",100,3000);
  gg1->Fit(fBILrat,"R");
  
  








  
  canvas->Print( (pdf_file_name + ")").c_str());

}

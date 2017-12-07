const int nmag = 5;
const double Ilinear[nmag] = {0.00307690156,0.0049948,0.00300464,0.004519309778,3.18840579710144957e-03}; //[GeV/A]
char const *magnets[nmag] = {"HB","Q1","Q2","Q3","Dipole"};

double fitResult(int kk, double P){

  if (kk==0){return -0.00330881+0.0014527*P+-7.72414e-05*pow(P,2)+3.83419e-06*pow(P,3)+2.47341e-05*pow(P,4)-3.36047e-06*pow(P,5)+1.38778e-07*pow(P,6);}
  else if (kk==1){return -0.000853647+0.0013583*P-0.000906733*pow(P,2)+0.000320065*pow(P,3)-5.81811e-05*pow(P,4)+4.77767e-06*pow(P,5)-1.17895e-07*pow(P,6);}
  else if (kk==2){return -3.64723e-07-6.68743e-08*P+2.96829e-07*pow(P,2)-1.02997e-07*pow(P,3)+1.29097e-08*pow(P,4)-6.04878e-10*pow(P,5)+5.83663e-12*pow(P,6);}
  else if (kk==3){return  -0.00110238+5.91188e-05*P-6.74454e-05*pow(P,2)+0.000189113*pow(P,3)-3.62886e-05*pow(P,4)+2.72258e-06*pow(P,5)-7.31832e-08*pow(P,6);}
  else {return -0.000878185-4.97453e-05*P+0.000481375*pow(P,2)-0.000183576*pow(P,3)+3.13612e-05*pow(P,4)+-2.54554e-06*pow(P,5)+7.96046e-08*pow(P,6);}
}


		 
void fitter(){
  //read input file
  TTree *tt = new TTree("tt","tt");
  tt->ReadFile("fitting.txt","pp:hb:q1:q2:q3:dip");
  const int npoints = tt->GetEntries();
  double t_P[npoints];
  double t_I[nmag][npoints];
  for (int jj=0;jj<npoints; jj++){
    tt->GetEntry(jj);
    TLeaf *t_PL = tt->GetLeaf("pp");
    TLeaf *t_IL_0 = tt->GetLeaf("hb");
    TLeaf *t_IL_1 = tt->GetLeaf("q1");
    TLeaf *t_IL_2 = tt->GetLeaf("q2");
    TLeaf *t_IL_3 = tt->GetLeaf("q3");
    TLeaf *t_IL_4 = tt->GetLeaf("dip");
    t_P[jj] = t_PL->GetValue();//GeV
    t_I[0][jj] = t_IL_0->GetValue();
    t_I[1][jj] = t_IL_1->GetValue();
    t_I[2][jj] = t_IL_2->GetValue();
    t_I[3][jj] = t_IL_3->GetValue();
    t_I[4][jj] = t_IL_4->GetValue();
  }

  //make output pdf
  TCanvas *canvas = new TCanvas("canvas","fit IvP SHMS", 800, 800);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin(0.2);
  //canvas->SetRightMargin(0.2);
  std::string pdf_file_name = "IvP_fits_shms.pdf";
  gROOT->SetBatch(true);
  canvas->Update();
  
  //setup the plots and fits
  TGraph *g1[nmag];
  TGraph *g2[nmag];
  TGraph *g3[nmag];
  TGraph *g4[nmag];
  TGraph *g5[nmag];

  TF1 *f1[nmag];
  TF1 *f3[nmag];
  double a[nmag][7];

  double linear[nmag][npoints];
  double linearDiff[nmag][npoints];
  double linearRes[nmag][npoints];
  double Iset_prime[nmag][npoints];
  double Iset_primeDiff[nmag][npoints];
  double Iset_primeRes[nmag][npoints];
  
  //loop magnets
  for (int ii=0; ii<nmag; ii++){
    if (ii==2) continue;

    // plot P vs I and fit
    g1[ii] = new TGraph(npoints,t_P,t_I[ii]);
    g1[ii]->SetMarkerStyle(22);
    g1[ii]->SetMarkerColor(kBlue);
    g1[ii]->SetMarkerSize(1.5);
    g1[ii]->SetTitle(Form("P vs I, %s",magnets[ii]));
    g1[ii]->GetXaxis()->SetTitle("P [GeV]");
    g1[ii]->GetYaxis()->SetTitle("Iset [A]");
    g1[ii]->GetYaxis()->SetTitleOffset(1.8);
    g1[ii]->Draw("ap");
    f1[ii] = new TF1(Form("f1_%d",ii),"pol1",t_P[0],t_P[npoints-1]);
    g1[ii]->Fit(f1[ii],"QR");
    canvas->Print( (pdf_file_name + "(").c_str());
    

    //plot Iset - Ilinear and fit
    for (int jj=0; jj<npoints; jj++){
      linear[ii][jj] = t_P[jj]/Ilinear[ii];
      linearDiff[ii][jj] = t_I[ii][jj]-linear[ii][jj];
      linearRes[ii][jj] = linearDiff[ii][jj]/linear[ii][jj];
    }
    
    g2[ii] = new TGraph(npoints,t_P,linearDiff[ii]);
    g2[ii]->SetMarkerStyle(22);
    g2[ii]->SetMarkerColor(kBlue);
    g2[ii]->SetMarkerSize(1.5);
    g2[ii]->SetTitle(Form("I_{set} - I_{linear}, %s",magnets[ii]));
    g2[ii]->GetXaxis()->SetTitle("P [GeV]");
    g2[ii]->GetYaxis()->SetTitle("I_{set} - I_{linear} [A]");
    g2[ii]->GetYaxis()->SetTitleOffset(1.8);
    g2[ii]->Draw("ap");
    canvas->Print( (pdf_file_name + "(").c_str());


    //plot rel residual (Iset-Ilinear)/Ilinear and fit
    g3[ii] = new TGraph(npoints,t_P,linearRes[ii]);
    g3[ii]->SetMarkerStyle(22);
    g3[ii]->SetMarkerColor(kBlue);
    g3[ii]->SetMarkerSize(1.5);
    g3[ii]->SetTitle(Form("I_{set} - I_{linear} / I_{linear}, %s",magnets[ii]));
    g3[ii]->GetXaxis()->SetTitle("P [GeV]");
    g3[ii]->GetYaxis()->SetTitle("I_{set} - I_{linear} / I_{set} [unitless]");
    g3[ii]->GetYaxis()->SetTitleOffset(1.8);
    g3[ii]->Draw("ap");
    f3[ii] = new TF1(Form("f3_%d",ii),"pol6",t_P[0],t_P[npoints-1]);
    g3[ii]->Fit(f3[ii],"QR");

    a[ii][0] = f3[ii]->GetParameter(0);
    a[ii][1] = f3[ii]->GetParameter(1);
    a[ii][2] = f3[ii]->GetParameter(2);
    a[ii][3] = f3[ii]->GetParameter(3);
    a[ii][4] = f3[ii]->GetParameter(4);
    a[ii][5] = f3[ii]->GetParameter(5);
    a[ii][6] = f3[ii]->GetParameter(6);

    cout<<"Magnet "<<ii<<":\t"<<a[ii][0]<<"+"<<a[ii][1]<<"P+"<<a[ii][2]<<"P^2+"<<a[ii][3]<<"P^3+"<<a[ii][4]<<"P^4+"<<a[ii][5]<<"P^5+"<<a[ii][6]<<"P^6"<<endl;
    
    canvas->Print( (pdf_file_name + "(").c_str());

    //plot Ilinear*(1+(Iset-Ilinear))/Ilinear = Iset'
    for (int jj=0; jj<npoints; jj++){
      //Iset_prime[ii][jj] = linear[ii][jj]*(1+g3[ii]->Eval(t_P[jj]));
      Iset_prime[ii][jj] = linear[ii][jj]*(1+fitResult(ii,t_P[jj]));
      
      Iset_primeDiff[ii][jj] = t_I[ii][jj]-Iset_prime[ii][jj];
      Iset_primeRes[ii][jj] = Iset_primeDiff[ii][jj]/t_I[ii][jj];
    }

    //plot Iset - Iset'
    g4[ii] = new TGraph(npoints,t_P,Iset_primeDiff[ii]);
    g4[ii]->SetMarkerStyle(22);
    g4[ii]->SetMarkerColor(kBlue);
    g4[ii]->SetMarkerSize(1.5);
    g4[ii]->SetTitle(Form("I_{set} - I'_{set} (from fit), %s",magnets[ii]));
    g4[ii]->GetXaxis()->SetTitle("P [GeV]");
    g4[ii]->GetYaxis()->SetTitle("I_{set} - I'_{set} [A]");
    g4[ii]->GetYaxis()->SetTitleOffset(1.8);
    g4[ii]->SetMinimum(-1);
    g4[ii]->SetMaximum(1);
    g4[ii]->Draw("ap");
    canvas->Print( (pdf_file_name + "(").c_str());

    //plot (Iset-Iset')/Iset (money plot, rel res)
    g5[ii] = new TGraph(npoints,t_P,Iset_primeRes[ii]);
    g5[ii]->SetMarkerStyle(22);
    g5[ii]->SetMarkerColor(kBlue);
    g5[ii]->SetMarkerSize(1.5);
    g5[ii]->SetTitle(Form("I_{set} - I'_{set} (from fit) / I_{set}, %s",magnets[ii]));
    g5[ii]->GetXaxis()->SetTitle("P [GeV]");
    g5[ii]->GetYaxis()->SetTitle("I_{set} - I'_{set} /I_{set} [unitless]");
    g5[ii]->GetYaxis()->SetTitleOffset(1.8);
    g5[ii]->SetMinimum(-0.0005);
    g5[ii]->SetMaximum(0.0005);
    g5[ii]->Draw("ap");
    canvas->Print( (pdf_file_name + "(").c_str());
    
  }//end loop mag

  canvas->Print( (pdf_file_name + ")").c_str());

}

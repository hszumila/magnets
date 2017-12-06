const double betaG = 0.009728667*0.999967; //units of [kG/A]
const double leffG = 1; //units of [m]
const double tol = 1E-5;
const double P0 = 0.004519309778; //units of [GeV/A]
const double increment = 1.0; //units of [A]

double beta(double Iset){

  //return betaG;

  //from root:
  return (5.57412E-24*pow(Iset,6) - 4.57913E-20*pow(Iset,5) + 1.34093E-16*pow(Iset,4) - 1.5001E-13*pow(Iset,3) + 1.28572E-12*pow(Iset,2) + 2.34223E-12*Iset + 0.00973899);

}

double leff(double Iset){

   return leffG;
  
}

double calcPratio(double I_linear, double I_iter){

  double ratio = (I_iter/I_linear)*(beta(I_iter)/betaG)*(leff(I_iter)/leffG);
  return ratio-1;

}


// input the desired momentum in [GeV]
void iteration(double P){

  //calculate the I_linear
  double Ilin = P/P0; //returns units of [A]

  double eta = 10000;
  int ii=0;
  double Iiter = Ilin;
  double bratio;
  double lratio;
  while (abs(eta)>tol){
    eta = calcPratio(Ilin, Iiter);
    //cout<<"\t"<<ii<<"\t"<<eta<<"\t"<<Iiter<<"\t"<<beta(Iiter)<<endl;
    if (eta > 0) {Iiter = Iiter*(1-abs(eta)/2.0);}
    else {Iiter = Iiter*(1+abs(eta)/2.0);}
    
    bratio = beta(Iiter)/betaG;
    lratio = leff(Iiter)/leffG;

    ii++;
    if (abs(ii) > 100) {break;}
  }


  cout<<"Requested energy:\t"<<P<<" [GeV]"<<endl;
  cout<<"Converged after "<<ii<<" iterations."<<endl;
  cout<<"Recommended Iset: "<<Iiter<<" [A], from initial guess: "<<Ilin<<" [A]."<<endl;
  cout<<"beta ratio: "<<bratio<<" leff ratio: "<<lratio<<endl;
  cout<<"Resulting difference from nominal: "<<eta<<endl;

}

const double betaG = 0.0118376*1.00088*0.999994; //units of [kG/A]
const double leffG = 1; //units of [m]
const double tol = 1E-5;
const double P0 = 3.18840579710144957e-03; //units of [GeV/A]
const double increment = 1.0; //units of [A]

double beta(double Iset){

  //return betaG;

  //from root:
  return (-9.72072E-25*pow(Iset,6) + 9.79429E-21*pow(Iset,5) - 3.80396E-17*pow(Iset,4) + 7.02145E-14*pow(Iset,3) - 5.80238E-11*pow(Iset,2) + 1.88857E-9*Iset + 0.0118584);

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

const double betaG = 0.0071009434*1.00035; //units of [kG/A]
const double leffG = 1.863; //units of [m]
const double tol = 1E-5;
const double P0 = 0.0049948; //units of [GeV/A]
const double increment = 1.0; //units of [A]

double beta(double Iset){

  //return betaG;

  //from root:
  return (2.24863E-23*pow(Iset,6) - 1.56067E-19*pow(Iset,5) + 3.72383E-16*pow(Iset,4) - 4.12478E-13*pow(Iset,3) + 2.29043E-10*pow(Iset,2) - 6.2683E-08*Iset + 0.00711028);

}

double leff(double Iset){

  // return leffG;

  return leffG* (1 + (Iset>1220)*(4.60576E-12*pow(Iset-1220,3)-1.00781E-8*pow(Iset-1220,2)+7.17211E-7*(Iset-1220)));

  //from dave:
  //return leffG*(1 -7.5934E-22*pow(Iset,6) + 7.9983E-18*pow(Iset,5) - 2.9479E-14*pow(Iset,4) + 4.7633E-11*pow(Iset,3) - 3.5855E-08*pow(Iset,2) + 1.1777E-05*Iset - 1.2738E-03);
  
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

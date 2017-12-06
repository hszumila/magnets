#include <stdlib.h>
#include "TMath.h"
#include "TROOT.h"

const double betaG = 8.94363583242782471E-03;//8.93508498183325471E-03*1.00096*0.999997; //units of [kG/A]
const double leffG = 5.28744767; //units of [m]
const double tol = 0.00001;//1E-5;
const double P0 = 3.26687082629409332E-03; //units of [GeV/A]
const double increment = 1.0; //units of [A]

double beta(double Iset){

  //from root:  
  //return (-8.58205E-23*pow(Iset,6) + 7.80459e-19*pow(Iset,5) - 2.57949e-15*pow(Iset,4) + 3.62541e-12*pow(Iset,3) - 2.11122e-09 *pow(Iset,2) + 3.00178e-07*Iset + 0.0090136 );

   return (7.76495E-29*pow(Iset,8) - 9.59365E-25*pow(Iset,7) +4.80895E-21*pow(Iset,6) -1.25206E-17*pow(Iset,5) +1.81193E-14*pow(Iset,4) - 1.48202E-11*pow(Iset,3) +6.77632E-9 *pow(Iset,2) -1.71809E-6*Iset + 0.00916971 );

}

double leff(double Iset){

  return leffG* (1 + (Iset>1220)*(2.5173E-12*pow(Iset-1220,3)-1.17767E-8*pow(Iset-1220,2)-3.60694E-6*(Iset-1220)));
  
}

double calcPratio(double I_linear, double I_iter){

  double ratio = (I_iter/I_linear)*(beta(I_iter)/betaG)*(leff(I_iter)/leffG);
  cout<<"ratio:\t"<<ratio<<endl;
  return ratio-1.0;

}


// input the desired momentum in [GeV]
void iteration(double P){

  //calculate the I_linear
  double Ilin = P/P0; //returns units of [A]

  double eta = 10000.0;
  int ii=0;
  double Iiter = Ilin;
  double Blin = betaG*Ilin;
  double Biter = Blin;
  double bratio;
  double lratio;
  
  while (fabs(eta)>tol){
    eta = calcPratio(Ilin, Iiter);
    cout<<"\t"<<ii<<"\t"<<eta<<"\t"<<Iiter<<"\t"<<beta(Iiter)/betaG<<"\ttol:\t"<<tol<<endl;
    if (eta > 0) {Iiter = Iiter*(1-fabs(eta)/2.0);}
    else {Iiter = Iiter*(1+fabs(eta)/2.0);}

    bratio = beta(Iiter)/betaG;
    lratio = leff(Iiter);

    Biter = (lratio * beta(Iiter) * Iiter)/10.0;

    ii++;
    if (ii > 100) {break;}   
  }
 
  
  cout<<"Requested energy:\t"<<P<<" [GeV]"<<endl;
  cout<<"Converged after "<<ii<<" iterations."<<endl;
  cout<<"Recommended NMR Bset: "<<Biter<<" [T], from initial guess: "<<Blin/10<<" [T]."<<endl;
  cout<<"Corresponding Iset: "<<Iiter<<" [A], from linear approx: "<<Ilin<<" [A]."<<endl;
  cout<<"beta ratio: "<<bratio<<" leff ratio: "<<lratio<<endl;
  cout<<"Resulting difference from nominal: "<<eta<<endl;

}

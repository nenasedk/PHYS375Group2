#include "Star.hh"

using namespace std;
Star(double dens, double temp, double aX, double aY, double aZ, double amu){
  central_dens = dens;
  central_temp = temp;
  _X = aX;
  _Y = aY;
  _Z = aZ;
  _mu = amu;

  R_0 = 1.0e-10;
  rk = new AdaptSolve();
}
// all functions must have a type (in most of our cases, void)
void Star::Mass(){ // Within the class, you can access all of the functions from wherever without needing to pass it as an argument
  //nvar is basically the order of the starting DE - specifically it is the number
  // of first order ODE's we must solve to fully solve our original DE. So it'll usually be 1 or 2
  int nvar = 1; //1 for mass because we only need to solve dm/dr (not d^2m/dr^2)
  vector<double> M(nvar); // the vector must have length
  
  // R_0 will always be the same - so we'll use a class variable R_0.
  // Likewise for the surface radius - we'll need to calculate that before we can evaluate some of these DEs.
  rk.RKSolve(M,nvar,R_0,R_surf,dMdr);

  // This is the routine for extracting the results from the integrator
  _Rad = vector<double>(rk.kount);// points along x where y was calculated
  _Rad= rk.xp; // This only needs to be done the FIRST TIME rk is used.
  
  _Massy = vector<vector<double> >(nvar,vector<double>(rk.kount)); // because we might be calculating multiple DEs, Massy is a vector of vectors, where each vecor is the solution for that nvar

  for(int i = 0; i<nvar;i++){// Get each of the solutions
    _Massy = rk.yp.at(i);
  }
  rk.reset()
    
}

void Star::Luminosity(){
  
  int nvar = 1;
  vector<double> L(nvar);
  rk.RKSolve(L,nvar,R_0,R_surf,dLdr);
  _Lum = L;
}

void Star::OptDepth(){
  int nvar = 1;
  vector<double> OptD(nvar);
  rk.RKSolve(OptD,nvar,R_0,R_surf,dtaudr);
  _OptD= OptD;
}

void Star::Temperature(){
  int nvar = 1;
  vector<double> Temp(nvar);

  rk = new AdaptSolve();
  rk.RKSolve(Temp,nvar,R_0,R_surf,dTdr);
  _Temp = Temp;
}

void Star::Density(){
  int nvar = 1;
  vector<double> Dens(nvar);
  double R_surf;   //x R_surf

  rk.RKSolve(Dens,nvar,R_0,R_surf,dpdr);
  _Dens = Dens;
}

// Partial Derivates of Pressure
//revised
void Star::dPdp(double R){//partial der of P wrt density

  double dP = (pow(3*pow(pi,2),2/3)/3)*(pow(Constants.hbar,2)/(Constants.me*Constants.mp))*pow(_Dens(R)/Constants.mp,2/3) + Constants.k*_Temp(R)/(mu*Constants.mp); // not sure how to call mu
  return dP;
}

//revised
void Star::dPdT(double R){
  
  double dP = _Dens(R)*Constants.k/(_mu*Constants.mp) + (4/3)*Constants.a*pow(_Temp(R),3);
  return dP;
  
}

//Energy Generation Rates 
//revised
void Star::EGR_PP(double R){ // will this function take the density and temperature as vectors or doubles?
  double dens_5 = _Dens(R)*1e-5;
  double T_6 = _Temp(R)*1e-6;
  double eps = 1.07e-7*dens_5*pow(_X,2)*pow(T_6,4); // not sure how to call X here
  return eps;
}
//revised
void Star::EGR_CNO(double R){// same as above fn but for CNO
  double dens_5 = _Dens(R)*1e-5;
  double T_6 = _Temp(R)*1e-6;
  double X_cno = 0.03*_X;
  double eps = 8.24e-26*dens_5*_X*X_cno*pow(T_6,19.9);
  return eps;
}

//Opacity Function
//revised
void Star::Opacity(double R){
  double dens_3 = _Dens(R)*1e-3;

  double Kes = 0.02*(1+_X);
  double Kff = 1.0e24*(_Z+0.0001)*pow(dens_3,0.7)*pow(_Temp(R),-3.5);
  double KH = 2.5e-32*(_Z/0.02)*pow(dens_3,0.5)*pow(_Temp(R),9);
  
  double OPsum = pow(KH,-1) + pow(max(Kes,Kff),-1);
  return pow(OPsum,-1);
}

//Pressure
//revised
void Star::Pressure(double R){
  
  double P = (pow(3*pow(pi,2),2/3)/5)*(pow(Constants.hbar,2)/(Constants.me)*pow(_Dens(R)/Constants.mp,5/3) + _Dens(R)*Constants.k*_Temp(R)/(_mu*Constants.mp) + (1/3)*Constants.a*pow(_Temp(R),4);
  return P;
				       }

//Derivatives wrt to r
// All derivative functions must have the form dYdX(double x, vector<double> y, vector<double> dydx)
// as arguments

void Star::dMdr(double R, vector<double>& M, vector<double>& dMdr){// mass change with radius
    dMdr = 4*pi*pow(R,2)*_Dens(R);// This also needs to be the density at radius R
    
}
//revised
void Star::dLdr(double R, vector<double>& L, vector<double>& dLdr){// luminosity change with radius

    dLdr = 4*pi*pow(R,2)*_Dens(R)*(Star.EGR_CNO(R) + Star.EGR_PP(R));

}
//revised
void Star::dtaudr(double R, vector<double>& OptD, vector<double>& dtaudr){

	dtaudr = Star.Opacity(R)*_Dens(R);

}
//revised
void Star::dTdr(double R, vector<double>& T, vector<double>& dTdr){
  
	dTdr = -min(3*Star.Opacity(R)*_Dens(R)*_Lum(R)/(16*pi*Constants.a*Constants.c*pow(_Temp(R),3)*pow(R,2)), (1-1/Constants.gamma)*_Temp(R)*Constants.G*_Mass(R)*_Dens(R)/(_Pres(R)*pow(R,2)));

}

//revised
void Star::dpdr(double R, vector<double>& Dens, vector<double> dpdr){// density change with radius

	dpdr = -(Constants.G*_Mass(R)*_Dens(R)*pow(R,-2) + Star.dPdT(R)*dTdr(R))/Star.dPdp(R); //unsure how to use dTdr here

}


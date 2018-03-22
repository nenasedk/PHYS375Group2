#include "Star.hh"

using namespace std;
Star(double dens, double temp, double aX, double aY, double aZ, double amu){
  central_dens = dens;
  central_temp = temp;
  X = aX;
  Y = aY;
  Z = aZ;
  mu = amu;

  R_0 = 1.0e-10;
  rk = new AdaptSolve();
}
// all functions must have a type (in most of our cases, void)
void Star::Mass(){ // Within the class, you can access all of the functions from wherever without needing to pass it as an argument

  //nvar is basically the order of the starting DE - specifically it is the number
  // of first order ODE's we must solve to fully solve our original DE. So it'll usually be 1 or 2
  int nvar = 1; //1 for mass because we only need to solve dm/dr (not d^2m/dr^2)
  vector<double> M(nvar); // the vector must have length

  double end;   //x end

  // R_0 will always be the same - so we'll use a class variable R_0.
  // Likewise for the surface radius - we'll need to calculate that before we can evaluate some of these DEs.
  rk.RKSolve(M,nvar,R_0,R_surf,dMdr);
  M = M;
}

Star::Luminosity(void dLdr){
  vector<double> L;
  int nvar;
  double end;   //x end
  rk.RKSolve(L,nvar,R_0,end,dLdr);
  M = M;
}

Star::OptDepth(void dtaudr){
  vector<double> OptD;
  int nvar;
  double end;   //x end

  rk.RKSolve(OptD,nvar,R_0,end,dtaudr);
  OptD= OptD;
}

Star::Temperature(void dTdr){
  vector<double> Temp;
  int nvar;
  double end;   //x end

  rk = new AdaptSolve();
  rk.RKSolve(Temp,nvar,R_0,end,dTdr);
  Temp = Temp;
}

Star::Density(void dpdr){
  vector<double> Dens;
  int nvar;
  double end;   //x end

  rk.RKSolve(Dens,nvar,R_0,end,dpdr);
  Dens = Dens;
}

// Partial Derivates of Pressure
Star::dPdp(double aDens, double aT){//partial der of P wrt density
  double dens = aDens
  double T = aT;
  double dP = (pow(3*pow(pi,2),2/3)/3)*(pow(Constants.hbar,2)/(Constants.me*Constants.mp)*pow(dens/Constants.mp,2/3) + Constants.k*T/(mu*Constants.mp); // not sure how to call mu
  return dP;
}

Star::dPdT(double aDens, double aT){
  double dens = aDens;
  double T = aT;
  
  double dP = dens*Constants.k/(mu*Constants.mp) + (4/3)*Constants.a*pow(T,3);
  return dP;
}

//Energy Generation Rates
Star::EGR_PP(double aDens, double aT){ // will this function take the density and temperature as vectors or doubles?
  double dens_5 = aDens*1e-5;
  double T_6 = aT*1e-6;
  double eps = 1.07e-7*dens_5*pow(X,2)*pow(T_6,4); // not sure how to call X here
  return eps;
}

Star::EGR_CNO(double aDens, double aT){// same as above fn but for CNO
  double dens_5 = aDens*1e-5;
  double T_6 = aT*1e-6;
  double X_cno = 0.03*X;
  double eps = 8.24e-26*dens_5*X*X_cno*pow(T_6,19.9);
  return eps;
}

//Opacity Function
Star::Opacity(double aX, double aZ, double aDens, double aT){
  double X = aX;
  double Z = aZ;
  double dens_3 = aDens*1e-3;
  double T = aT;
  
  double Kes = 0.02*(1+X);
  double Kff = 1.0e24*(Z+0.0001)*pow(dens_3,0.7)*pow(T,-3.5);
  double KH = 2.5e-32*(Z/0.02)*pow(dens_3,0.5)*pow(T,9);
  
  double OPsum = pow(KH,-1) + pow(max(Kes,Kff),-1);
  return pow(OPsum,-1);
}

//Pressure
Star::Pressure(double aDens, double aT){
  double dens = aDens;
  double T = aT;
  
  double P = (pow(3*pow(pi,2),2/3)/5)*(pow(Constants.hbar,2)/(Constants.me)*pow(dens/Constants.mp,5/3) + dens*Constants.k*T/(mu*Constants.mp) + (1/3)*Constants.a*pow(T,4);
  return P;
				       }

//Derivatives wrt to r
// All derivative functions must have the form dYdX(double x, vector<double> y, vector<double> dydx)
// as arguments

    void Star::dMdr(double R, vector<double> M, vector<double> dMdr){// mass change with radius
    dMdr = 4*pi*pow(R,2)*dens;// This also needs to be the density at radius R
    
  }

Star::dLdr(double ar, double aDens , double aT, vector<double> dydx){// luminosity change with radius
  double r = ar;
  double dens = aDens;
  double T = aT;
  vector<double> dL = 4*pi*pow(r,2))*dens*(Star.EGR_CNO(dens, T) + Star.EGR_PP(dens, T));
  return dL;
}

Star::dtaudr(double ar, double aDens , double aT, vector<double> dydx){
  double r = ar;
  double dens = aDens;
  double T = aT;
  vector<double> dtau = Star.Opacity(X,Z,dens,T)*dens;
  return dtau;
}

Star::dTdr(double ar, double aDens, double aL, double aT, double aM, double aP){
  double r = ar;
  double dens = aDens;
  double L = aL;
  double T = aT;
  double M = aM;
  double P = aP;
  
  vector<double> dT = -min(3*Star.Opacity(X,Z,dens,T)*dens*L/(16*pi*Constants.a*Constants.c*pow(T,3)*pow(r,2)), (1-1/Constants.gamma)*T*Constants.G*M*dens/(P*pow(r,2)));
  return dT;
}

Star::dpdr(double ar, double aDens, double aM, double aT, vector<double> dydx){// density change with radius
  double r = ar;
  double dens = aDens;
  double M = aM;
  double T = aT;
  vector<double> dp = -(Constants.G*M*dens*pow(r,-2) +Star.dPdT(dens, T)*Star.dTdr())/Star.dPdp(dens, T);
  return dp;
}


#include "Star.hh"

using namespace std;

Star::Star(double dens, double temp, double aX, double aY, double aZ, double amu){
  central_dens = dens;
  central_temp = temp;
  _X = aX;
  _Y = aY;
  _Z = aZ;
  _mu = amu;

  R_0 = 1.0e-10;
  rk = AdaptSolve();
}

void Star::Reset(){
  fill(_Dens.begin(), _Dens.end(), 0.0);
  fill(_Temp.begin(), _Temp.end(), 0.0);
  fill(_Mass.begin(), _Mass.end(), 0.0);
  fill(_Lum.begin(), _Lum.end(), 0.0);
  fill(_OptD.begin(), _OptD.end(), 0.0);
  fill(_Rad.begin(), _Rad.end(), 0.0);
  fill(_Pres.begin(), _Pres.end(), 0.0);
  rk.Reset();
}

// y = state vector
// y = {Density, Temperature, Mass, Luminosity, Optical depth}
void Star::Derivatives(double x, vector<double> &y, vector<double> &dydx){
  // Density
  dydx.at(0) = dpdr(x,y.at(0),y.at(2),y.at(1),dydx.at(1)); //need to fix units
  
  // Temperature
  dydx.at(1) = dTdr(x,y.at(0),y.at(2),y.at(1),y.at(3));

  // Mass
  dydx.at(2) = dMdr(x,y.at(1));
  
  // Luminosity
  dydx.at(3) = dLdr(x,y.at(0),y.at(1));

  // Optical depth
  dydx.at(4) = dtaudr(y.at(1),y.at(2));
  
}

void Star::EvaluateAll(){
  vector<double> state = vector<double>(5,0.0);
  int nvar = 5;
  rk.RKSolve(state,nvar,R_0,R_surf,Derivatives);
  _Rad = vector<double>(rk.kount);
  _Rad = rk.xp;

  _Dens = rk.yp.at(0);
  _Temp = rk.yp.at(1);
  _Mass = rk.yp.at(2);
  _Lum = rk.yp.at(3);
  _OptD = rk.yp.at(4);
  
  int len = rk.kount;
  for(int i = 0; i<len;i++){
    _Pres.at(i) = Pressure(_Rad.at(i),_Dens.at(i),_Temp.at(i));
  }
}
    

double Star::dPdp(double aDens, double aT){//partial der of P wrt density
  double dP = (pow(3*pow(M_PI,2),2/3)/3)*(pow(hbar,2))/(me*mp)*pow(aDens/mp,2/3) + k*aT/(_mu*mp); // not sure how to call mu
  return dP;
}

double Star::dPdT(double R,double dens,double temp){
  double dP = dens*k/(_mu*mp) + (4/3)*a*pow(temp,3);
  return dP;
  
}



  
//Energy Generation Rates 
//revised
double Star::EGR_PP(double R, double dens, double temp){ // will this function take the density and temperature as vectors or doubles?
  double dens_5 = dens*1e-5;
  double T_6 = temp*1e-6;
  double eps = 1.07e-7*dens_5*pow(_X,2)*pow(T_6,4); // not sure how to call X here
  return eps;
}
//revised
double Star::EGR_CNO(double R, double dens, double temp){// same as above fn but for CNO
  double dens_5 = dens*1e-5;
  double T_6 = temp*1e-6;
  double X_cno = 0.03*_X;
  double eps = 8.24e-26*dens_5*_X*X_cno*pow(T_6,19.9);
  return eps;
}
double Star::EGR_3a(double R, double dens, double temp){// same as above fn but for CNO
  double dens_5 = dens*1e-5;
  double T_8 = temp*1e-8;
  double eps = 3.85e-8*pow(dens_5,2)*pow(_Y,3)*pow(T_8,44.0);
  return eps;
}



//Opacity Function
//revised
double Star::Opacity(double dens, double temp){
  double dens_3 = dens*1e-3;

  double Kes = 0.02*(1+_X);
  double Kff = 1.0e24*(_Z+0.0001)*pow(dens_3,0.7)*pow(temp,-3.5);
  double KH = 2.5e-32*(_Z/0.02)*pow(dens_3,0.5)*pow(temp,9);
  
  double OPsum = pow(KH,-1) + pow(max(Kes,Kff),-1);
  return pow(OPsum,-1);
}

//Pressure
//revised
double Star::Pressure(double R,double dens,double temp){
  
  double P = (pow(3*pow(M_PI,2),2/3)/5)*(pow(hbar,2)/(me)*pow(dens/mp,5/3) +
				       dens*k*temp/(_mu*mp) + (1/3)*a*pow(temp,4));
  return P;
}


//Derivatives wrt to r
// Mass change with radius
// This also needs to be the density at radius R
double Star::dMdr(double R,double dens){
  double dM = 4*M_PI*pow(R,2)*dens;
  return dM;
}
//revised
double Star::dLdr(double R, double dens, double temp){// luminosity change with radius
  double dL = 4*M_PI*pow(R,2)*dens*(EGR_CNO(R,dens,temp) + EGR_PP(R,dens,temp) + EGR_3a(R,dens,temp));
  return dL;
}
//revised
double Star::dtaudr(double dens, double temp){
  double dtau = Opacity(dens,temp)*dens;
  return dtau;
}
//revised
double Star::dTdr(double R, double dens, double mass, double temp, double lum){
  double rad  = 3.0 *Opacity(dens,temp)*dens*lum / (16*M_PI*a*c*pow(temp,3)*pow(R,2));
  double conv = (1. - 1.0/agamma)* temp*G*mass*dens/(Pressure(R,dens,temp)*pow(R,2));
  return min(rad,conv);
}
    
// density change with radius
double Star::dpdr(double R,double dens, double temp, double mass,double dt){
  double dp = -(G*mass*dens*pow(R,-2) + dPdT(R,dens,temp)*dt)/dPdp(dens,temp); 
  return dp;

}


#include "Star.hh"

using namespace std;
Star(double dens, double temp, double aX, double aY, double aZ, double amu){
  central_dens = dens;
  central_temp = temp;
  X = aX;
  Y = aY;
  Z = aZ;
  mu = amu;
}

Star::y(){
  vector<double> y;
  int nvar;

  double start; //x(0)
  double end;   //x end

  rk = new AdaptSolve();
  rk.RKSolve(y,nvar,start,end,dydx);
  pres = y;
}

Star::dydx(double x, vector<double>y , vector<double> dydx){
  
}


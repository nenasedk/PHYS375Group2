#include "Star.hh"
#include "AdaptSolve.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "math.h"
using namespace std;
void EvaluateAll(Star*,AdaptSolve *);
void Derivatives(double, vector<double>&,vector<double>&);
double Dens = 0.0;
double Temp = 0.0;
double X = 0.0;
double Y = 0.0;
double Z = 0.0;
double mu = 0.0;

Star *s = new Star(Dens,Temp,X,Y,Z,mu);
void EvaluateAll(Star *s, AdaptSolve *rk ){

  vector<double> state = vector<double>(5,0.0);
  state.at(0) = s->central_dens; // Density
  state.at(1) = s->central_temp; // Temperature
  state.at(2) = 0.0; // Mass
  state.at(3) = 0.0; // Luminosity
  state.at(4) = s->Opacity(s->central_dens,s->central_temp)*s->central_dens; // Opacity
  int nvar = 5;
  
  rk->RKSolve(state,nvar,s->R_0,10.0,&Derivatives);
  s->_Rad = rk->xp;
  s->_Dens = rk->yp.at(0);
  s->_Temp = rk->yp.at(1);
  s->_Mass = rk->yp.at(2);
  s->_Lum = rk->yp.at(3);
  s->_OptD = rk->yp.at(4);
  int len = rk->xp.size();
  vector<double> pres = vector<double>(len);
  for(int i = 0; i<len;i++){
    pres.at(i) = s->Pressure(s->_Rad.at(i),s->_Dens.at(i),s->_Temp.at(i));
  }
  s->_Pres = pres;
  cout << "Evaluated a star!" << endl;
  
}

// y = state vector
// y = {Density, Temperature, Mass, Luminosity, Optical depth}
void Derivatives(double x, vector<double> &y, vector<double> &dydx){
  // Density
  dydx.at(0) = s->dpdr(x,y.at(0),y.at(1),y.at(2),dydx.at(1)); //need to fix units
  
  // Temperature
  dydx.at(1) = s->dTdr(x,y.at(0),y.at(1),y.at(2),y.at(3));

  // Mass
  dydx.at(2) = s->dMdr(x,y.at(0));
  
  // Luminosity
  dydx.at(3) = s->dLdr(x,y.at(0),y.at(1));

  // Optical depth
  dydx.at(4) = s->dtaudr(y.at(0),y.at(1));
  
}


int main(){
  AdaptSolve *rk = new AdaptSolve();
  for(int loop = 1; loop < 101; loop++){
    double Temp = 1.0e6;
    double Dens = 1.5e5;
    X = 0.734;
    Y = 0.250;
    Z = 0.016;
    mu = pow((2.0*X + 0.75*Y + 0.5*Z),-1);
    s->NewStar(Dens, Temp, X, Y, Z, mu);
    EvaluateAll(s,rk);
    ostringstream fileName;
    fileName << "MSStar_" << loop << ".txt";
	
    //cout << fileNameResult << endl;
    ofstream myfile (fileName.str().c_str());
    cout << "Writing Star " << loop << " to file." << endl;
    if (myfile.is_open()){
      for (int i = 0; i < s->_Rad.size(); i++){
	myfile << s->_Rad.at(i) << "," << s->_Dens.at(i) << "," << s->_Temp.at(i) << "," << s->_Mass.at(i) << "," << s->_Lum.at(i) << "," << s->_OptD.at(i) << endl;	
      }
    }

    myfile.close();
    rk->Reset();
  }
  delete rk;
  return 0;
}

#include "Star.hh"
#include "AdaptSolve.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "math.h"
using namespace std;
void EvaluateAll(Star*,AdaptSolve *,double,double,int,int,double);
void Derivatives(double, vector<double>&,vector<double>&);

// Gross global variables that we have to use 
double Dens = 1.0e6;
double Temp = 1.5e5;
double X = 0.7;
double Y = 0.28;
double Z = 0.02;
double mu = 2.4;

// Because I can't figure out a better way to pass the deriv function
// - can't be a class member, but needs access to star functions
Star *s = new Star(Dens,Temp,X,Y,Z,mu);




// EvaluateAll
// This function takes a star class and a de solver, and it fills all of our state variables
void EvaluateAll(Star *aStar, AdaptSolve *rk,double Rad,double h,int nsave,int maxstep,double dx){
  // The state vector:
  vector<double> state = vector<double>(6,0.0);
  state.at(0) = aStar->central_dens; // Density
  state.at(1) = aStar->central_temp; // Temperature
  state.at(2) = 1e-10;             // Mass
  state.at(3) = 1e-10;             // Luminosity
  state.at(4) = aStar->Opacity(aStar->central_dens,aStar->central_temp)*aStar->central_dens; // Optical depth
  state.at(5) = 100.0; // Opacity proxy for BCs
  int nvar = 6; // Size of state vector

  
  aStar->R_surf = Rad;
  // Set up the de solver
  rk->SetStep(h); // Initial step size to try
  rk->SetNSave(nsave); // How many data points we want saved
  rk->SetMaxSteps(maxstep); // Maximum number of steps  to integrate (10 million is usually safe)
  rk->SetSaveInterval(dx); // How far apart we want our R values saved
  rk->RKSolve(state,nvar,aStar->R_0,aStar->R_surf,&Derivatives); //BCs are taken care of in the solver

  // Store variables. Turns out to actually be unnecessary, but oh well.
  aStar->_Rad = rk->xp;
  aStar->_Dens = rk->yp.at(0);
  aStar->_Temp = rk->yp.at(1);
  aStar->_Mass = rk->yp.at(2);
  aStar->_Lum = rk->yp.at(3);
  aStar->_OptD = rk->yp.at(4);
  int len = rk->xp.size();
  vector<double> pres = vector<double>(len);
  for(int i = 0; i<len;i++){
    pres.at(i) = aStar->Pressure(aStar->_Rad.at(i),aStar->_Dens.at(i),aStar->_Temp.at(i));
  }
  aStar->_Pres = pres;
  //cout << "Evaluated a star!" << endl;
  
}

/* Derivative function
 * Calls of of the equations of stellar structure from the star class,
 * and has the correct arguments for the de solver
 *
 * Note, do not try to make this a class member function. (because pointers)
 * y = state vector
 * y = {Density, Temperature, Mass, Luminosity, Optical depth}
 */
void Derivatives(double x, vector<double> &y, vector<double> &dydx){
  // Density
  dydx.at(0) = s->dpdr(x,y.at(0),y.at(1),y.at(2),dydx.at(1));  
  // Temperature
  dydx.at(1) = s->dTdr(x,y.at(0),y.at(1),y.at(2),y.at(3));
  // Mass
  dydx.at(2) = s->dMdr(x,y.at(0));  
  // Luminosity
  dydx.at(3) = s->dLdr(x,y.at(0),y.at(1));
  // Optical depth
  dydx.at(4) = s->dtaudr(y.at(0),y.at(1));
  // BCs
  dydx.at(5) = s->OpBC(y.at(0),y.at(1),dydx.at(1));  
}

// Bisection method for finding central density
Star* Bisection(Star *a, Star *b, Star *c){
  double eps = 1e-2;
  int i = 0; // limit number of trials
  AdaptSolve *as = new AdaptSolve();
  cout << "Bisection method to tune density..." << endl;
  while((b->central_dens - a->central_dens)/2.0 > eps && i<30){
    if(a->LumBisec() * c->LumBisec() < 0.0){
      b = c;
    }
    else{
      a = c;
    }
    double middens = (b->central_dens + a->central_dens)/2.0;
    as->SetConvergence(1.0);
    c->NewStar(middens,c->central_temp,c->_X,c->_Y,c->_Z,c->_mu);
    EvaluateAll(c,as,1.0e10,1.0e4,1e5,1e6,5.0e4);
    as->Reset();
    i++;
  }
  EvaluateAll(c,as,1.0e10,1.0e4,1e5,1e6,5.0e4);
  delete as;
  return c;  
}

// Run the program
int main(){
  AdaptSolve *rk = new AdaptSolve();
  Star *a = new Star(Dens,Temp,X,Y,Z,mu);
  Star *b = new Star(Dens,Temp,X,Y,Z,mu);
  for(int loop = 1; loop < 101; loop++){
    // Initial Conditions
    Temp = 0.1 * loop * 1.0e6; //Linearly scaling the central temperature
    Dens = 2.0e5;
    X = 0.734;
    Y = 0.250;
    Z = 0.016;
    mu = pow((2.0*X + 0.75*Y + 0.5*Z),-1);
    a->NewStar(0.9*Dens,Temp,X,Y,Z,mu);
    b->NewStar(1.1*Dens,Temp,X,Y,Z,mu);
    // Set up our star and evaluate
    s->NewStar(Dens, Temp, X, Y, Z, mu);
    EvaluateAll(a,rk,1.0e10,1.0e4,10000,10000000,5.0e4);
    rk->Reset();
    EvaluateAll(b,rk,1.0e10,1.0e4,10000,10000000,5.0e4);
    rk->Reset();
    EvaluateAll(s,rk,1.0e10,1.0e4,10000,10000000,5.0e4);
    s = Bisection(a,b,s);
    cout << "Evaluated a star!" << endl;
    // File output
    ostringstream fileName;
    fileName << "Data/MSStar_" << loop << ".txt";

    ofstream myfile (fileName.str().c_str());
    cout << "Writing Star " << loop << " to file." << endl;
    if (myfile.is_open()){
      int i = 0;
      // File header
      myfile << "X = " << X << endl;
      myfile << "Y = " << Y << endl;
      myfile << "Z = " << Z << endl;
      myfile << "mu = " << mu << endl;
      myfile << "Radius" << "," << "Density" << "," << "Temp"  << "," << "Mass" << "," << "Lum" << "," <<"OptD" << "," << "Pres"<< endl;
      // Data out
      while (rk->xp.at(i)<rk->xp.at(i+1)){
	myfile << rk->xp.at(i) << "," << rk->yp.at(0).at(i) << "," << s->_Temp.at(i) << "," << s->_Mass.at(i) << "," << s->_Lum.at(i) << "," << s->_OptD.at(i) << "," << s->_Pres.at(i)<< endl;
	//myfile << rk->xp.at(i) << "," <<  rk->yp.at(2).at(1) << endl;
	i++;
	
      }
    }

    myfile.close();
    rk->Reset();
  }
  delete rk;
  delete a;
  delete b;
  delete s;
  return 0;
}

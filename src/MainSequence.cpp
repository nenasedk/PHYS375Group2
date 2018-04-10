#include "Star.hh"
#include "AdaptSolve.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "math.h"
#include <iomanip>
using namespace std;
Star EvaluateAll(Star,AdaptSolve *,double,double,int,int,double);
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
Star s(Dens,Temp,X,Y,Z,mu);

// EvaluateAll
// This function takes a star class and a de solver, and it fills all of our state variables
Star EvaluateAll(Star aStar, AdaptSolve *rk,double Rad,double h,int nsave,int maxstep,double dx){
  // The state vector:
  vector<double> state = vector<double>(6,0.0);
  state.at(0) = aStar.central_dens; // Density
  state.at(1) = aStar.central_temp; // Temperature
  state.at(2) = 4./3. * M_PI * pow(aStar.R_0,3.0)*aStar.central_dens;             // Mass
  state.at(3) = aStar.dLdr(aStar.R_0,aStar.central_dens,aStar.central_temp)*aStar.R_0/3.0; // Luminosity
  state.at(4) = aStar.dtaudr(aStar.central_dens,aStar.central_temp)*aStar.R_0; // Optical depth
  state.at(5) = 100.0; // Opacity proxy for BCs -IC doesn't matter just must be >1
  int nvar = 6; // Size of state vector
  aStar.R_surf = Rad;
  
  // Set up the de solver
  rk->SetStep(h); // Initial step size to try
  rk->SetNSave(nsave); // How many data points we want saved
  rk->SetMaxSteps(maxstep); // Maximum number of steps  to integrate (10 million is usually safe)
  rk->SetSaveInterval(dx); // How far apart we want our R values saved
  s = Star(aStar);
  int err = rk->RKSolve(state,nvar,0.01,s.R_surf,&Derivatives); //BCs are taken care of in the solver
  // Store variables. Turns out to actually be unnecessary, but oh well.
  aStar._Rad = rk->xp;
  
  aStar._Dens = rk->yp.at(0);
  aStar._Temp = rk->yp.at(1);
  aStar._Mass = rk->yp.at(2);
  aStar._Lum = rk->yp.at(3);
  aStar._OptD = rk->yp.at(4);
  aStar._dLdr = s._dLdr;
  aStar._PP = s._PP;
  aStar._CNO = s._CNO;
  aStar._3a = s._3a;
  s.Reset();
  return aStar;  
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
  // Luminosity
  dydx.at(3) = s.dLdr(x,y.at(0),y.at(1));  
  // Temperature
  dydx.at(1) = s.dTdr(x,y.at(0),y.at(1),y.at(2),y.at(3));  
  // Density
  dydx.at(0) = s.dpdr(x,y.at(0),y.at(1),y.at(2),dydx.at(1));
  // Mass
  dydx.at(2) = s.dMdr(x,y.at(0));  
  // Optical depth
  dydx.at(4) = s.dtaudr(y.at(0),y.at(1));
  // BCs
  dydx.at(5) = s.OpBC(y.at(0),y.at(1),dydx.at(0));  
}

/* Bisection method for finding central density
 *
 * The bisection method is a root-finding method
 * that iteratively shrinks the bounds on either side
 * of the zero, until the difference is within 
 * a specified tolerance
 */
Star Bisection(Star aStar, Star bStar, Star cStar){
  double eps = 1.e-2; // Allowed tolerance on variation in density
  int i = 0; // limit number of trials (30)
  // Check to ensure that the 0 exists within our bounds
  if(aStar.LumBisec() * bStar.LumBisec() > 0.0){
    cout << "No Root within bounds, exiting." << endl;
    cStar.Reset();
    return cStar;
  }
  
  AdaptSolve *as = new AdaptSolve();
  cout << "Bisection method to tune density..." << endl;
  while(fabs(bStar.central_dens - aStar.central_dens)/2.0 > eps){
    double al  = aStar.LumBisec();
    double bl = bStar.LumBisec();
    double cl = cStar.LumBisec();
    if(i == 30){
      // Check if we've run too long, and pick a non-radiative star
      if(aStar._Temp.at(aStar.SurfRad()) < aStar._Temp.at(aStar.SurfRad())){
	cStar.NewStar(aStar.central_dens,cStar.central_temp,cStar._X,cStar._Y,cStar._Z,cStar._mu);
      } else {
	cStar.NewStar(bStar.central_dens,cStar.central_temp,cStar._X,cStar._Y,cStar._Z,cStar._mu);
      }
      cout << "Choosing lower star, exiting bisection" << endl;
      break;
    }

    // Bisection check for luminosity error
    if(al * cl < 0.0){
      bStar = Star(cStar);
    }
    else{
      aStar = Star(cStar);
    }

    // Find the new average to test
    double middens = (bStar.central_dens + aStar.central_dens)/2.0;
    cStar.central_dens = middens;
    cStar = EvaluateAll(cStar,as,1.0e10,5.0e4,10000,10000000,5.0e4);
    i++;
    as->Reset();
  }

  // Ensure that we pick a star with a sensible surface temp
  if(aStar._Temp.at(aStar.SurfRad()) > 5.0e5){
    cStar.central_dens = bStar.central_dens;
  } else if (bStar._Temp.at(aStar.SurfRad()) > 5.0e5){
    cStar.central_dens = aStar.central_dens;
  }
  cout << "Found a star!" <<endl;
  cStar = EvaluateAll(cStar,as,1.0e10,5.0e4,10000,10000000,5.0e4);
  cStar.SurfRad();
  cStar.FillOpacity();
  cStar.FillPres();
  cStar.FillEGR();
  cStar.FilldPdT();
  cStar.FilldLdr();
  cStar.FilldTdr();
  aStar.Reset();
  bStar.Reset();
  delete as;
  return cStar;  
}

// Run the program
int main(){
  AdaptSolve *rk = new AdaptSolve();
  Star a(Dens,Temp,X,Y,Z,mu);
  Star b(Dens,Temp,X,Y,Z,mu);
  Star c(Dens,Temp,X,Y,Z,mu);
  int nstar = 100;
  for(int loop = 13; loop < 20; loop++){
    // Initial Conditions
    Temp = pow(10.,((7.5-6.2)/(float(nstar))*loop + 6.2)); //Power Law scaling the central temperature
    X = 0.734;
    Y = 0.250;
    Z = 0.016;
    mu = pow((2.0*X + 0.75*Y + 0.5*Z),-1);
    
    a.NewStar(300.,Temp,X,Y,Z,mu);
    b.NewStar(5.0e5,Temp,X,Y,Z,mu); 
    c.NewStar((5.003e5/2.0), Temp, X, Y, Z, mu);
    
    // Set up our star and evaluate
    a = EvaluateAll(a,rk,1.0e10,5.0e4,10000,10000000,5.0e4);
    rk->Reset();
    b = EvaluateAll(b,rk,1.0e10,5.0e4,10000,10000000,5.0e4);
    rk->Reset();
    c = EvaluateAll(c,rk,1.0e10,5.0e4,10000,10000000,5.0e4);
    c = Bisection(a,b,c);
    cout << "Evaluated a star!" << endl;
    // File output
    if(c._Dens.size() == 0){
      continue;
    }
    ostringstream fileName;
    fileName << "DataNewTemps7/MSStar_" << loop << ".txt";

    ofstream myfile (fileName.str().c_str());
    cout << "Writing Star " << loop << " to file." << endl;
    if (myfile.is_open()){
      int i = 0;
      // File header
      myfile << "X = " << X << endl;
      myfile << "Y = " << Y << endl;
      myfile << "Z = " << Z << endl;
      myfile << "mu = " << mu << endl;
      myfile << setw(10) << "Radius" << ", " << setw(10) << "Density" << ", " << setw(10) << "Temp"  << ", " << setw(10) << "Mass" << ", " << setw(10) << "Lum" << ", " << setw(10) << "OptD" << ", " << setw(10) << "Kff"  << ", " << setw(10) << "KH"  << ", " << setw(10) << "Kes"  << ", " << setw(10) << "PP"  << ", " << setw(10) << "CNO"  << ", " << setw(10) << "3a"  << ", " << setw(10) << "Deg"  << ", " << setw(10) << "Gas"  << ", " << setw(10) << "Rad" << ", " << setw(10) << "TotPres"<< ", " << setw(10) << "dLdr" << ", " << setw(10) << "dPdT, " << setw(10) << "dTRad" << ", " << setw(10) << "dTConv" << endl;
      // Data out
      while (i<=c._MaxRad){
	myfile << setw(10) << c._Rad.at(i) << ", ";
	myfile << setw(10) << c._Dens.at(i) << ", ";
	myfile << setw(10) << c._Temp.at(i) << ", ";
	myfile << setw(10) << c._Mass.at(i) << ", ";
	myfile << setw(10) << c._Lum.at(i) << ", " ;
	myfile << setw(10) << c._OptD.at(i) << ", ";
	myfile << setw(10) << c._Kff.at(i) << ", ";
	myfile << setw(10) << c._KH.at(i) << ", ";
	myfile << setw(10) << c._Kes.at(i) << ", ";
	myfile << setw(10) << c._PP.at(i) << ", ";
	myfile << setw(10) << c._CNO.at(i) << ", ";
	myfile << setw(10) << c._3a.at(i) << ", ";
	myfile << setw(10) << c._DegPres.at(i) << ", ";
	myfile << setw(10) << c._GasPres.at(i) << ", ";
	myfile << setw(10) << c._RadPres.at(i) << ", ";
	myfile << setw(10) << c._Pres.at(i) << ", ";
	myfile << setw(10) << c._dLdr.at(i) << ", ";
	myfile << setw(10) << c._dPdT.at(i) << ", ";
	myfile << setw(10) << c._dTRad.at(i) << ", ";
	myfile << setw(10) << c._dTConv.at(i) << endl;
	//myfile << rk->xp.at(i) << "," <<  rk->yp.at(2).at(1) << endl;
	i++;
	
      }
    }

    myfile.close();
    a.Reset();
    b.Reset();
    c.Reset();
    rk->Reset();
  }

  cout << "You did it! Good Job! :)"<< endl;
  delete rk;
  return 0;
}

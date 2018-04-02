#include "Star.hh"
#include "AdaptSolve.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "math.h"
using namespace std;
Star EvaluateAll(Star*,AdaptSolve *,double,double,int,int,double);
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
  state.at(3) = aStar.dLdr(aStar.R_0,aStar.central_dens,aStar.central_temp)*aStar.R_0/3.0;             // Luminosity
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
  //cout << s.central_dens << ", " <<  s.R_0 << endl;
  int err = rk->RKSolve(state,nvar,s.R_0,s.R_surf,&Derivatives); //BCs are taken care of in the solver
  // Store variables. Turns out to actually be unnecessary, but oh well.
  aStar._Rad = rk->xp;
  
  aStar._Dens = rk->yp.at(0);
  aStar._Temp = rk->yp.at(1);
  aStar._Mass = rk->yp.at(2);
  aStar._Lum = rk->yp.at(3);
  aStar._OptD = rk->yp.at(4);
  /*aStar._DegPres = s._DegPres;
  aStar._RadPres = s._RadPres;
  aStar._GasPres = s._GasPres;
  aStar._Kff = s._Kff;
  aStar._Kes = s._Kes;
  aStar._KH = s._KH;*/
  aStar._PP = s._PP;
  aStar._CNO = s._CNO;
  aStar._3a = s._3a;
  //aStar._Pres = s._Pres;
  s.Reset();
  return aStar;
  return aStar;
  
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
  dydx.at(0) = s.dpdr(x,y.at(0),y.at(1),y.at(2),dydx.at(1));  
  // Temperature
  dydx.at(1) = s.dTdr(x,y.at(0),y.at(1),y.at(2),y.at(3));
  // Mass
  dydx.at(2) = s.dMdr(x,y.at(0));  
  // Luminosity
  dydx.at(3) = s.dLdr(x,y.at(0),y.at(1));
  // Optical depth
  dydx.at(4) = s.dtaudr(y.at(0),y.at(1));
  // BCs
  dydx.at(5) = s.OpBC(y.at(0),y.at(1),dydx.at(0));  
}

// Bisection method for finding central density
Star Bisection(Star aStar, Star bStar, Star cStar){
  double eps = 1.e-3;
  int i = 0; // limit number of trials
  AdaptSolve *as = new AdaptSolve();
  cout << "Bisection method to tune density..." << endl;
  while(fabs(bStar.central_dens - aStar.central_dens)/2.0 > eps && i<50){
    //cout << "Test" << endl;
    //cout << i << endl;
    if(aStar.LumBisec() * cStar.LumBisec() < 0.0){
      //cout << "test 1" << endl;
      //b.NewStar(c.central_dens,c.central_temp,c._X,c._Y,c._Z,c._mu);
      //b = EvaluateAll(b,as,1.0e10,1.0e4,10000,10000000,5.0e4);
      bStar = Star(cStar);
    }
    else{
      //cout << "test 2" << endl;
      //a.NewStar(c.central_dens,c.central_temp,c._X,c._Y,c._Z,c._mu);
      //a = EvaluateAll(a,as,1.0e10,1.0e4,10000,10000000,5.0e4);
      aStar = Star(cStar);
    }
    //cout << aStar._OptD.at(3000) << ", " << bStar._OptD.at(3000) << ", " << cStar._OptD.at(3000) << endl;
    as->Reset();
    double middens = (bStar.central_dens + aStar.central_dens)/2.0;
    //as.SetConvergence(0.001);
    //cout << middens << endl;
    cStar.central_dens = middens;
    //Star temp(middens,cStar.central_temp,cStar._X,cStar._Y,cStar._Z,cStar._mu);
    cStar = EvaluateAll(cStar,as,1.0e10,1.0e3,1000,10000000,1.0e5);
    //cStar = Star(temp);
    //cout << cStar._Rad.at(3000) << endl;
    //as->Reset();
    i++;
    if(i==48){cout << "Bisection method did not converge" << endl;}
  }
  cout << "Found a star!" <<endl;
  //as->SetConvergence(1.0e-4);
  cStar = EvaluateAll(cStar,as,1.0e10,1.0e3,10000,10000000,1.0e5);
  cStar.FillOpacity();
  cStar.FillPres();
  delete as;
  return cStar;  
}

// Run the program
int main(){
  AdaptSolve *rk = new AdaptSolve();
  Star a(Dens,Temp,X,Y,Z,mu);
  Star b(Dens,Temp,X,Y,Z,mu);
  Star c(Dens,Temp,X,Y,Z,mu);
  for(int loop = 34; loop < 35; loop++){
    // Initial Conditions
    //Temp = 2.0e5*loop + 5.0e6; //Linearly scaling the central temperature
    Temp = pow(10.,0.00909091*loop +6.6); //Power Law scaling the central temperature
    Dens = 1.1e5;
    X = 0.734;
    Y = 0.250;
    Z = 0.016;
    mu = pow((2.0*X + 0.75*Y + 0.5*Z),-1);
    
    a.NewStar(0.7*Dens,Temp,X,Y,Z,mu);
    b.NewStar(20.0*Dens,Temp,X,Y,Z,mu); 
    c.NewStar(Dens, Temp, X, Y, Z, mu);
    
    // Set up our star and evaluate
    a = EvaluateAll(a,rk,1.0e10,1.0e3,1000,10000000,5.0e5);
    //cout << rk->yp.at(4).size() << ", " << a._OptD.size() << endl;
    rk->Reset();
    b = EvaluateAll(b,rk,1.0e10,1.0e3,1000,10000000,5.0e5);
    rk->Reset();
    c = EvaluateAll(c,rk,1.0e10,1.0e3,1000,10000000,5.0e5);
    c = Bisection(a,b,c);
    cout << "Evaluated a star!" << endl;
    // File output
    ostringstream fileName;
    fileName << "DataNewTemps5/MSStar_" << loop << ".txt";

    ofstream myfile (fileName.str().c_str());
    cout << "Writing Star " << loop << " to file." << endl;
    if (myfile.is_open()){
      int i = 0;
      // File header
      myfile << "X = " << X << endl;
      myfile << "Y = " << Y << endl;
      myfile << "Z = " << Z << endl;
      myfile << "mu = " << mu << endl;
      myfile << "Radius" << "," << "Density" << "," << "Temp"  << "," << "Mass" << "," << "Lum" << "," <<"OptD" << "," << "Kff"  << "," << "KH"  << "," << "Kes"  << "," << "PP"  << "," << "CNO"  << "," << "3a"  << "," << "Deg"  << "," << "Gas"  << "," << "Rad" << "," << "TotPres"<< endl;
      // Data out
      while (c._Rad.at(i) < c._Rad.at(i+1)){
	myfile << c._Rad.at(i) << ",";
	myfile << c._Dens.at(i) << ",";
	myfile << c._Temp.at(i) << ","  ;
	myfile << c._Mass.at(i) << ",";
	myfile << c._Lum.at(i) << ",";
	myfile << c._OptD.at(i) << ",";
	myfile << c._Kff.at(i) << ",";
	myfile << c._KH.at(i) << ",";
	myfile << c._Kes.at(i) << ",";
	myfile << c._PP.at(i) << ",";
	myfile << c._CNO.at(i) << ",";
	myfile << c._3a.at(i) << ",";
	myfile << c._DegPres.at(i) << ",";
	myfile << c._GasPres.at(i) << ",";
	myfile << c._RadPres.at(i) << ",";
	myfile << c._Pres.at(i)<< endl;
	//myfile << rk->xp.at(i) << "," <<  rk->yp.at(2).at(1) << endl;
	i++;
	
      }
    }

    myfile.close();
    rk->Reset();
  }
  cout << "Closing..."<< endl;
  delete rk;
  return 0;
}

/*
 * Star Class
 *
 * This class will define the PDEs used to describe stellar structure and 
 * and evolution, and will allow metallicity to be defined as an 
 * independant variable to allow for variation of that parameter
 *
 * It will call the AdaptMesh class to solve the PDEs, and will output 
 * the data necessary to make the required plots (plots will be done in python)
 */

#ifndef _Star_
#define _Star_

#include <cmath>
#include <iostream>
#include <vector>
#include "Constants.hh"
#include "AdaptSolve.hh"
#include <algorithm>
using namespace std;

class Star{

public:
  // densc,Tc,X,Y,Z,mu
  Star(long double,long double,long double,long double,long double,long double);
  ~Star();
  //void EvaluateAll();
  //void Derivatives(long double,vector<long double>&,vector<long double>&);
  void Reset();
  void NewStar(long double,long double,long double,long double,long double,long double);

  // I should really make get functions for all of these...
  vector<long double> _Rad; // vector of evaluation points, defined by SetSaveInterval in AdaptSolve
  vector<long double> _Mass; // mass at location x
  vector<long double> _Temp; // temperature at location x
  vector<long double> _Dens; // density at location x
  vector<long double> _Lum; // luminosity at location x
  vector<long double> _OptD; // optical depth at location x
  vector<long double> _Pres; // Pressure at x
  
  long double central_dens; 
  long double central_temp;
  // I'm using _'s just as a convention for denoting class level variables
  long double _X,_Y,_Z;
  long double _mu;
  
  //AdaptSolve *rk;

  long double R_0, R_surf; // starting radius, surface radius

  long double Pressure(long double,long double,long double);
  long double EGR_PP(long double,long double,long double);
  long double EGR_CNO(long double,long double,long double);
  long double EGR_3a(long double,long double,long double);
  long double Opacity(long double,long double);

  long double dPdp(long double,long double);
  long double dPdT(long double,long double,long double);
  long double OpBC(long double,long double,long double);
  long double dpdr(long double,long double,long double,long double,long double);
  long double dTdr(long double,long double,long double,long double,long double);
  long double dMdr(long double, long double);
  long double dLdr(long double, long double, long double);
  long double dtaudr(long double, long double);
  long double LumBisec();
  int SurfRad();
  int MaxArg();
private:
  void init();

};

#endif

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
#include <stdexcept>
using namespace std;

class Star{

public:
  // densc,Tc,X,Y,Z,mu
  Star(double,double,double,double,double,double);
  Star(const Star &);
  ~Star();
  //void EvaluateAll();
  //void Derivatives(double,vector<double>&,vector<double>&);
  void Reset();
  void NewStar(double,double,double,double,double,double);

  // I should really make get functions for all of these...
  vector<double> _Rad; // vector of evaluation points, defined by SetSaveInterval in AdaptSolve
  vector<double> _Mass; // mass at location x
  vector<double> _Temp; // temperature at location x
  vector<double> _Dens; // density at location x
  vector<double> _Lum; // luminosity at location x
  vector<double> _OptD; // optical depth at location x
  vector<double> _Pres; // Pressure at x
  vector<double> _dLdr;
  vector<double> _dPdT;
  
  vector<double> _DegPres; // Pressure at x
  vector<double> _GasPres;
  vector<double> _RadPres; // Pressure at x
  
  vector<double> _Kff; // luminosity at location x
  vector<double> _KH; // optical depth at location x
  vector<double> _Kes; // Pressure at x

  vector<double> _PP; // luminosity at location x
  vector<double> _CNO; // optical depth at location x
  vector<double> _3a; // Pressure at x
  
  double central_dens; 
  double central_temp;
  // I'm using _'s just as a convention for denoting class level variables
  double _X,_Y,_Z;
  double _mu;
  int _MaxRad;
  
  //AdaptSolve *rk;

  double R_0, R_surf; // starting radius, surface radius

  double Pressure(double,double,double);
  double EGR_PP(double,double,double);
  double EGR_CNO(double,double,double);
  double EGR_3a(double,double,double);
  double Opacity(double,double);

  double dPdp(double,double);
  double dPdT(double,double,double);
  double OpBC(double,double,double);
  double dpdr(double,double,double,double,double);
  double dTdr(double,double,double,double,double);
  double dMdr(double, double);
  double dLdr(double, double, double);
  double dtaudr(double, double);
  double LumBisec();
  int SurfRad();
  int MaxArg();
  void FillPres();
  void FillOpacity();
  void FillEGR();
  void FilldPdT();
  void FilldLdr();
private:
  void init();

};

#endif

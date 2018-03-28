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

using namespace std;
class Star{

public:
  // densc,Tc,X,Y,Z,mu
  Star(double,double,double,double,double);
  ~Star();
  void init();
  void Pressure(double);
  void EGR_PP(double);
  void EGR_CNO(double);
  void Opacity(double);
  void Density();
  void Temperature();
  void Luminosity();
  void Mass(void);
  void OptDepth();
  void dPdp(double);
  void dPdT(double);
  void dpdr(double, vector<double>,vector<double>);
  void dTdr(double, vector<double>,vector<double>);
  void dMdr(double, vector<double>,vector<double>);
  void dLdr(double, vector<double>,vector<double>);
  void dtaudr(double, vector<double>,vector<double>);
  
  void Loop();
  void Output();
private:
  double central_dens; 
  double central_temp;
  // I'm using _'s just as a convention for denoting class level variables
  double _X,_Y,_Z;
  double _mu;
  vector<double> _Pres;
  AdaptSolve::AdaptSolve rk;

  double R_0, R_surf; // starting radius, surface radius
  vector<double> _Rad; // vector of evaluation points, defined by SetSaveInterval in AdaptSolve
  vector<vector<double> > _Mass; // mass at location x
  vector<vector<double> > _Temp; // temperature at location x
  vector<vector<double> > _Dens; // density at location x
  vector<vector<double> > _Lum; // luminosity at location x
  vector<vector<double> > _OptD; // optical depth at location x
};

#endif

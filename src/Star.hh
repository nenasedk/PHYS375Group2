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
  void Pressure(double, double);
  void EGR_PP(double, double);
  void EGR_CNO(double, double);
  void Density();
  void Temperature();
  void Luminosity();
  void Mass(void);
  void OptDepth();
  void Opacity(double, double, double, double);
  void dPdp(double, double);
  void dPdT(double, double);
  void dpdr();
  void dTdr();
  void dMdr(double, vector<double>,vector<double>);
  void dLdr();
  void dtaudr();
  
  void Loop();
  void Output();
private:
  double central_dens; 
  double central_temp;
  double X,Y,Z;
  double mu;
  vector<double> pres;
  AdaptSolve::AdaptSolve rk;

  double R_0, R_surf; // starting radius, surface radius
  double _Mass;
};

#endif

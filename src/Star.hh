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
  void Pressure();
  void dPdr(double, vector<double>&, vector<double>&);
  
  void Loop();
  void Output();
private:
  double central_dens,central_temp,X,Y,Z,mu;
  vector<double> pres;
};

#endif

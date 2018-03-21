/*
 * Adaptive Mesh Class
 *
 * This class will define the adaptive mesh DE solver used to solve the Stars
 * PDEs. An implementation of RK4 will also be included. 
 *
 * This will provide the numerical component of the Stars project
 */

#ifndef _AdaptSolve_
#define _AdaptSolve_

#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

class AdaptSolve{

 public:
  AdaptSolve();
  ~AdaptSolve();
  void init();
  void SetConvergence(double);
  void RKSolve(vector<double>&, int, double, double, 
	       void (*derivs)(double, vector<double>&, vector<double>&));
  void SetSaveInterval(double);
  void Reset();
 private:
  
  double eps;
  double h1;
  double hmin;
  int *nok;
  int *nbad;
  int kount;
  int kmax;
  double step;
  int f_nvar;
  double f_h;
  // Outputs
  vector<double> xp;
  vector<vector<double> > yp;
  double dxsav;
  /* Butcher Tableau (Wikipedia)
   * To specify a particular method, one needs to provide the integer s 
   * (the number of stages), and the coefficients bij (for 1 ≤ j < i ≤ s), 
   * ci (for i = 1, 2, ..., s) and ai (for i = 2, 3, ..., s). The matrix [bij] 
   * is called the Runge–Kutta matrix, while the bi and ci are known as the 
   * weights and the nodes. 
   */
  static const double
    a2 = 1.0/5.0,  b21 = 1.0/5.0,
    a3 = 3.0/10.0, b31 = 3.0/40,     b32 = 9.0/40.0,
    a4 = 6.0/10.0, b41 = 3.0/10,     b42 = -9.0/10.0, b43 = 1.2,
    a5 = 1.0,      b51 = -11.0/54.0, b52 = 2.5,       b53 = -70.0/27.0, b54 = 35.0/27.0,
    a6 = 0.875,    b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0,
    b64 = 44275.0/110592.0, b65 = 253.0/4096.0,

    c1 = 37.0/278.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    dc5 = -277.0/14336.0;

  double dc1, dc3, dc4, dc6;
 
  void rkqs( vector<double>&,  vector<double>&, int, double *, double,
	     vector<double>&, double *, double *,
	    void (*)(double,  vector<double>&,  vector<double>&));

  void rkck( vector<double>&,  vector<double>&, int, double, double,  vector<double>&,
	     vector<double>&, void (*)(double,  vector<double>&,  vector<double>&));
 
};

#endif

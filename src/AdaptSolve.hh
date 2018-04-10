/*
 * Adaptive Mesh Class
 *
 * This class defines the adaptive mesh DE solver used to solve the Stars
 * PDEs. It is based of the RCKC method from Numerical Recipes in C. 
 *
 * This will provide the numerical component of the Stars project
 */

#ifndef _AdaptSolve_
#define _AdaptSolve_

#include <math.h>
#include <vector>
#include <iostream>
#include <stdexcept>
using namespace std;

//using namespace std:
// Class declaration
class AdaptSolve{
public:
  // Adapt Solve
  // Class Constructor - Creates an instance of the class when called
  AdaptSolve();
  
  // ~AdaptSolve
  // Class Destructor - deletes and closes class
  ~AdaptSolve();

  // Set functions
  void SetConvergence(double);
  void SetStep(double);
  void SetNSave(int);
  void SetMaxSteps(int);
  // RKSolve
  // Based on Numerical Recipes in C Chapter 16
  // Fourth-Order, Adaptive step size runge kutta integrator
  //
  // This is the actual solver function that we will use for integration
  // It takes in a vector address y, an int number of variables, double start, double end,
  //               and a function that is the derivative of y, which in turn takes arguments of
  //               double x, vector address of y, vector address of dydx.

  int RKSolve(vector<double>&, int, double, double, 
	       void (*)(double, vector<double>&, vector<double>&));

  // Set Save Interval
  void SetSaveInterval(double);

  // Reset
  // Clear all class variables so we can start again
  void Reset();

  // Outputs
  vector<double> xp; // Saved x values
  vector<vector<double> > yp; // Saved y values
  int kount;   //


  
  // Private variables and internal functions
private: 
  double eps;  // convergence value
  double hmin; // Minimum step size
  double hmax;
  int nok;    // Number of good integration steps
  int nbad;   // Number of bad integration steps (not converged)
  int kmax;    //
  double step; //
  int f_nvar;  // Number of variables
  double f_h;  // Current step size

  double dxsav; // How often to save
  int f_maxstep;
  
  /* Butcher Tableau (Wikipedia)
   * To specify a particular method, one needs to provide the integer s 
   * (the number of stages), and the coefficients bij (for 1 ≤ j < i ≤ s), 
   * ci (for i = 1, 2, ..., s) and ai (for i = 2, 3, ..., s). The matrix [bij] 
   * is called the Runge–Kutta matrix, while the bi and ci are known as the 
   * weights and the nodes. 
   */
  static const double
  a2 = 1.0/5.0,  b21 = 1.0/5.0,
    a3 = 3.0/10.0, b31 = 3.0/40.0,     b32 = 9.0/40.0,
    a4 = 6.0/10.0, b41 = 3.0/10.0,     b42 = -9.0/10.0,   b43 = 1.2,
    a5 = 1.0,      b51 = -11.0/54.0,   b52 = 2.5,         b53 = -70.0/27.0, b54 = 35.0/27.0,
    a6 = 0.875,    b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0,
    b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
    
    c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    dc5 = 277.0/14336.0;
  
  double dc1, dc3, dc4, dc6;

  // Initialises variables to sensible defaults
  void init();

  bool BCs(double,vector<double>&,vector<double>&);
  // Quality controlled Runge Kutta - makes sure global errors don't accumulate
  void rkqs( vector<double>&,  vector<double>&, int, double *, double,
	     vector<double>&, double*, double*,
	     void (*)(double,  vector<double>&,  vector<double>&));
  
  // Runge Kutta Cash - Karl Method
  // Fifth order in error so we can do adaptive step size to fourth order
  void rkck( vector<double>&,  vector<double>&, int, double, double,  vector<double>&,
	     vector<double>&, void (*)(double,  vector<double>&,  vector<double>&));
};
#endif

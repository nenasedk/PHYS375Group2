#include "AdaptSolve.hh"
// We only need to include this .hh file because the rest are included in the .hh file
using namespace std;

// Class Constructors
AdaptSolve::AdaptSolve(){
  init();
}
AdaptSolve::~AdaptSolve(){}

// Initialisation
void AdaptSolve::init(){
  eps = 0.0001;
  f_h = 0.01;
  hmin = 1.0e-9;
  hmax = 5000000.;
  kmax = 100000;
  f_maxstep = 1000000;

  // Dependant Butcher table variables
  dc1=    c1- (2825.0/27648.0);
  dc3=    c3- (18575.0/48384.0);
  dc4=    c4- (13525.0/55296.0);
  dc6=    c6-0.25;
  dxsav = 10.0;
}

// Set and reset functions
void AdaptSolve::SetStep(double h){f_h = h;}
void AdaptSolve::SetNSave(int s){kmax = s;}
void AdaptSolve::SetMaxSteps(int m){f_maxstep = m;}
void AdaptSolve::SetConvergence(double ep){eps = ep;}
void AdaptSolve::SetSaveInterval(double k){dxsav = k;}
void AdaptSolve::Reset(){
  kount = 0;
  nok = 0;
  nbad = 0;
  f_nvar = 5;
  f_h = 0.01;
  xp.clear();
  yp.clear();  
}


/* Adapt Solve
 *
 * This is the adaptive RK4 solver for a generic set of ODEs.
 *
 * Based on the examples given in Numerical Recipes in C, The Art of Scientific Computing 2ed, CH16
 *
 * ystart: the initial values of each of the DEs to be solved
 * nvar: the number of DEs to be solved
 * x1: starting position
 * x2: final position
 * derivs: a function that finds dydx for a given x and y.
 * 
 */
int AdaptSolve::RKSolve(vector<double>& ystart, int nvar, double x1, double x2, 
			 void (derivs)(double,vector<double>&, vector<double>&)){
  int maxstp = f_maxstep; // How many steps to take
  double tiny = 1e-30; // Prevent underflow
  int a = nvar; // How many DEs are we solving
  int i = 0;
  double x = 0.0;

  // Step sizes
  double hnext = 0.1;
  double hdid = 0.1;
  double h = f_h;
  double xsav = 0.0;
  
  vector<double> yscal = vector<double>(nvar,0.0);
  vector<double> y = vector<double>(nvar,0.0);
  vector<double> dydx = vector<double>(nvar,0.0);

  xp = vector<double>(maxstp);
  yp = vector<vector<double> >(nvar,vector<double>(maxstp,0.0));

  // Initialise starting position, step size, and ICs
  x=x1;
  h = (x2-x1) >= 0.0 ? fabs(f_h) : -fabs(f_h);
  nok = (nbad) = kount = 0;
  bool check = true;
  for (i=0;i<nvar;i++) {
    y.at(i)=ystart.at(i);
  }
  if (kmax > 0) xsav=x-dxsav*2.0; 

  // Main loop to fill xp and yp from x1 to x2
  for (int nstp=0; nstp<maxstp; nstp++) { //Take at most maxstp steps.
    (*derivs)(x,y,dydx);
    // Scaling used to monitor accuracy.
    for (i=0;i<nvar;i++){
      yscal.at(i)=fabs(y.at(i))+fabs(dydx.at(i)*h)+tiny;
    }
    //Store data.
    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
      xp.at(++kount) = x; 
      for (i=0;i<nvar;i++) yp.at(i).at(kount) = y.at(i);
      xsav=x;
    }
    
    // Check that we're not going to overshoot, and if so, it's probably a fully radiative star
    if (x+h - x2 > 0.0){
      cout << "Purely radiative star, reached integration limit" <<endl;
      return -1;}

    // Call the Adaptive Step size function
    rkqs(y,dydx,nvar,&x,h,yscal,&hdid,&hnext,derivs);
    // Check if this step was successful
    if (hdid == h) ++(nok); else ++(nbad);
    // Check completion
    if (BCs(x,y,dydx)) { //Are we done? 
      if (kmax) {
	xp.at(kount) = x; 
	for (i=0;i<nvar;i++) yp.at(i).at(kount) = y.at(i);
      }
      return 0; 
    }
    h = hnext;

    // Check that our step size isn't too small
    if (fabs(hnext) <= hmin){
      throw out_of_range("Step size too small in AdaptSolve");
      //h= hmin;
    }
    h = hnext;
  }
  cout << "Too many steps in routine AdaptSolve"<<endl;
}

/* Boundary Conditions
 *
 * This function returns true if any of the boundary conditions are satisfied
 *
 * Current BCs include:
 * Optical depth proxy
 * Density below 0 check
 * High mass check
 * Large radius check
 */
bool AdaptSolve::BCs(double x, vector<double>& y,vector<double>& dydx){
  if(y.at(0)<0.0){
    cout << "End on Dens BC" << endl;
    return true;}
  if(dydx.at(5)< 1.e-2){
    return true;}// dydx 5 is the opacity BC
  if(y.at(2) > 2e33){
    cout << "End on Mass BC" << endl;
    return true;}
  if(x>1.0e12){return true;}
  return false;
}

/* Quality step RK solver
 *
 * This function uses the RKCK method, and adapts the current step size
 * iteratively until the error on this step is within the specified tolerance.
 *
 * If the current error is smaller than the tolerance, the step size is allowed to grow.
 */
void AdaptSolve::rkqs(vector<double>& y, vector<double>& dydx, int n, double *x, double htry,
		      vector<double>& yscal, double *hdid, double *hnext,
		      void (*derivs)(double, vector<double>&, vector<double>&)){
  // Constants for safe step size change
  double safe = 0.9;
  double grow = -0.5;
  double shrink = -0.25;
  double errcon = 1.89e-2; 
  double errmax,h,htemp,xnew;
  
  vector<double> yerr= vector<double>(n);
  vector<double> ytemp= vector<double>(n);
  h=htry;

  // Loop until within specified error tolerance eps
  for (;;) {
    
    // Call Cash-Karl Runge Kutta
    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
    // Check if within tolerance
    errmax=0.0;
    for (int i=0;i<n-2;i++){
      double diff = fabs(yerr.at(i));
      errmax=max(errmax,diff/fabs(yscal.at(i)));
    }
    
    // Check completion
    if (errmax/eps <= 1.0) {break;}

    // Shrink step size
    h = h *  min(max(safe * pow((eps/(errmax + 1e-30)),-1.*shrink), 0.5 ),2.0);
    xnew=(*x)+h;
    if (xnew == *x){
      throw out_of_range("stepsize underflow in rkqs");
    } 
  }

  // Grow step size
  *hnext = h * min( max(safe * pow((eps/(errmax + 1e-30)),-1.*grow), 0.5 ), 2. );

  // Save and output
  if(*hnext > hmax){*hnext=hmax;}
  *x += (*hdid=h);
  for (int i=0;i<n;i++){ y.at(i)=ytemp.at(i);}
}

/* Runge Kutta Cash Karp method
 *
 * This function implements the fourth order RCKC method
 * with 5th order error correction
 *
 * This method takes in a state vector of length n,
 * and solves the set of n coupled or uncoupled DEs
 * given by the user provided derivs function
 * at position x
 */
void AdaptSolve::rkck(vector<double>& y, vector<double>& dydx, int n, double x,
		      double h, vector<double>& yout,
		      vector<double>& yerr, void (*derivs)(double, vector<double>&, vector<double>&)){
  double errmax,htemp,xnew;
  vector<double>ak2 = vector<double>(n);
  vector<double>ak3 = vector<double>(n);
  vector<double>ak4 = vector<double>(n);
  vector<double>ak5 = vector<double>(n);
  vector<double>ak6 = vector<double>(n);
  vector<double> ytemp = vector<double>(n);
  // 0th
  for(int i=0; i<n; i++){
    ytemp.at(i) = y.at(i) + b21*h*dydx.at(i);
  }
  // 1st
  for (int i=0;i<n;i++){
    ytemp.at(i) = y.at(i) + b21*h*dydx.at(i);
  }
  if(ytemp.at(0) < 0.0){
    ytemp.at(0) = y.at(0) + h*dydx.at(0);
  }
  (*derivs)(x+a2*h,ytemp,ak2);
  //2nd
  for (int i=0;i<n;i++){
    ytemp.at(i) = y.at(i) + h*(b31*dydx.at(i) + b32*ak2.at(i));
  }
  if(ytemp.at(0) < 0.0){
    ytemp.at(0) = y.at(0) + h*dydx.at(0);
  }
  (*derivs)(x+a3*h,ytemp,ak3);
  //3rd
  for (int i=0;i<n;i++){
    ytemp.at(i) = y.at(i) + h*(b41*dydx.at(i) + b42*ak2.at(i) + b43*ak3.at(i));
  }
  if(ytemp.at(0) < 0.0){
    ytemp.at(0) = y.at(0) + h*dydx.at(0);
  }
  (*derivs)(x+a4*h,ytemp,ak4);
  // 4th
  for (int i=0;i<n;i++){
    ytemp.at(i) = y.at(i) + h*(b51*dydx.at(i) + b52*ak2.at(i) + b53*ak3.at(i) + b54*ak4.at(i));
  }
  if(ytemp.at(0) < 0.0){
    ytemp.at(0) = y.at(0) + h*dydx.at(0);
  }
  (*derivs)(x+a5*h,ytemp,ak5);
  // 5th
  for (int i=0;i<n;i++){
    ytemp.at(i) = y.at(i) + h*(b61*dydx.at(i) + b62*ak2.at(i) + b63*ak3.at(i) + b64*ak4.at(i) + b65*ak5.at(i));
  }
  if(ytemp.at(0) < 0.0){
    ytemp.at(0) = y.at(0) + h*dydx.at(0);
  }
  (*derivs)(x+a6*h,ytemp,ak6);

  // Output y and error
  for (int i=0;i<n;i++){
    yout.at(i) = y.at(i) + h*(c1*dydx.at(i) + c3*ak3.at(i) + c4*ak4.at(i) + c6*ak6.at(i));
  }
  for (int i=0;i<n;i++){
    yerr.at(i) = h*(dc1*dydx.at(i) + dc3*ak3.at(i) + dc4*ak4.at(i) + dc5*ak5.at(i) + dc6*ak6.at(i));
  }
}

#include "AdaptSolve.hh"

using namespace std;

AdaptSolve::AdaptSolve(){
  init();
}
AdaptSolve::~AdaptSolve(){}

void AdaptSolve::init(){
  eps = 1e-8;
  h1 = 0.1;
  hmin = 1e-31;
  *nok = 0;
  *nbad = 0;
  kmax = 64;
  
  dc1=c1-2825.0/27648.0;
  dc3=c3-18575.0/48384.0;
  dc4=c4-13525.0/55296.0;
  dc6=c6-0.25;
  dxsav = 0.00001;
}

void AdaptSolve::SetConvergence(double ep){eps = ep;}
void AdaptSolve::SetSaveInterval(double k){dxsav = k;}
void AdaptSolve::Reset(){
  fill(xp.begin(), xp.end(), 0);
  fill(yp.begin(), yp.end(), 0);
}


/* Adapt Solve
 *
 * This is the adaptive RK4 solver for a generic set of ODEs.
 *
 * Based on the examples given in Numerical Recipes in C, The Art of Scientific Computing 2ed, CH16
 * 
 */
void AdaptSolve::RKSolve(vector<double>& ystart, int nvar, double x1, double x2, 
			  void (*derivs)(double,vector<double>&, vector<double>&)){
  int maxstp = 100000;
  double tiny = 1e-30;

  int a = nvar;
  
  int nstp = 0;
  int i = 0;
  double xsav = 0.0;
  double x = 0.0;
  double hnext = 0.0;
  double hdid = 0.0;
  double h = 0.0;
  
  vector<double> yscal= vector<double>(nvar);
  vector<double> y= vector<double>(nvar);
  vector<double> dydx= vector<double>(nvar);

  xp = vector<double>(nvar);
  yp = vector<vector<double> >(nvar,vector<double>(nvar));
  
  x=x1;
  h = (x2-x1) >= 0.0 ? fabs(h1) : -fabs(h1);
  *nok = (*nbad) = kount = 0;


  for (i=1;i<=nvar;i++) { y.at(i)=ystart.at(i); }
  if (kmax > 0) xsav=x-dxsav*2.0; 
  for (nstp=1;nstp<=maxstp;nstp++) { //Take at most maxstp steps.
    (*derivs)(x,y,dydx);
    for (i=1;i<=nvar;i++)
      //Scaling used to monitor accuracy. This general-purpose choice can be modified
      //if need be.
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+tiny;

    
    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
      xp[++kount]=x; //Store intermediate results.
      for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
      xsav=x;
    }

    
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;

    // Call RKQS
    rkqs(y,dydx,nvar,&x,h,yscal,&hdid,&hnext,derivs);
    // ***
    if (hdid == h) ++(*nok); else ++(*nbad);
    if ((x-x2)*(x2-x1) >= 0.0) { //Are we done?
      for (i=1;i<=nvar;i++){ ystart.at(i)=y.at(i);}
      if (kmax) {
	xp[++kount]=x; //Save final step.
	for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
      }
      delete &dydx;
      delete &y;
      delete &yscal;
      return; //Normal exit.
    }
    if (fabs(hnext) <= hmin) cout <<"Step size too small in AdaptSolve"<<endl;
    h=hnext;
  }
  cout << "Too many steps in routine AdaptSolve"<<endl;
}






void AdaptSolve::rkqs(vector<double>& y, vector<double>& dydx, int n, double *x, double htry,
		      vector<double>& yscal, double *hdid, double *hnext,
		      void (*derivs)(double, vector<double>&, vector<double>&)){
  double safe = 0.9;
  double grow = -0.2;
  double shrink = -0.25;
  double errcon = 1.89e-4; 
  double errmax,h,htemp,xnew;
  
  vector<double> yerr= vector<double>(n);
  vector<double> ytemp= vector<double>(n);
  h=htry; 
  for (;;) {
    // Call Cash-Karl Runge Kutta
    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
    
    errmax=0.0;  
    for (int i=1;i<=n;i++){
      errmax=max(errmax,fabs(yerr[i]/yscal.at(i)));
    }
    errmax /= eps; 
    if (errmax <= 1.0) {break;}
    htemp=safe*h*pow(errmax,shrink);
    h=(h >= 0.0 ? max(htemp,0.1*h) : min(htemp,0.1*h));
    xnew=(*x)+h;
    if (xnew == *x){
      cout << "stepsize underflow in rkqs"<< endl;
    }
  }
  
  if (errmax > errcon){
    *hnext=safe*h*pow(errmax,grow);}
  
  else{*hnext=5.0*h;}
  
  *x += (*hdid=h);
  for (int i=1;i<=n;i++){ y.at(i)=ytemp[i];}
  
  delete &ytemp;
  delete &yerr;
}


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

  for(int i = 1; i<=n; i++){
    ytemp.at(i) = y.at(i) + b21*h*dydx.at(i);
  }
  for (int i=1;i<=n;i++)
    ytemp.at(i)=y.at(i)+b21*h*dydx.at(i);
  (*derivs)(x+a2*h,ytemp,ak2); 
  for (int i=1;i<=n;i++)
    ytemp.at(i)=y.at(i)+h*(b31*dydx.at(i)+b32*ak2.at(i));
  (*derivs)(x+a3*h,ytemp,ak3);
  for (int i=1;i<=n;i++)
    ytemp.at(i)=y.at(i)+h*(b41*dydx.at(i)+b42*ak2.at(i)+b43*ak3.at(i));
  (*derivs)(x+a4*h,ytemp,ak4);
  for (int i=1;i<=n;i++)
    ytemp.at(i)=y.at(i)+h*(b51*dydx.at(i)+b52*ak2.at(i)+b53*ak3.at(i)+b54*ak4.at(i));
  (*derivs)(x+a5*h,ytemp,ak5); 
  for (int i=1;i<=n;i++)
    ytemp.at(i)=y.at(i)+h*(b61*dydx.at(i)+b62*ak2.at(i)+b63*ak3.at(i)+b64*ak4.at(i)+b65*ak5.at(i));
  (*derivs)(x+a6*h,ytemp,ak6); 
  for (int i=1;i<=n;i++)
    yout.at(i)=y.at(i)+h*(c1*dydx.at(i)+c3*ak3.at(i)+c4*ak4.at(i)+c6*ak6.at(i));
  for (int i=1;i<=n;i++)
    yerr.at(i)=h*(dc1*dydx.at(i)+dc3*ak3.at(i)+dc4*ak4.at(i)+dc5*ak5.at(i)+dc6*ak6.at(i));
}
  
    

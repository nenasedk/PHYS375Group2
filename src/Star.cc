#include "Star.hh"

using namespace std;

Star::Star(double dens, double temp, double aX, double aY, double aZ, double amu){
  central_dens = dens;
  central_temp = temp;
  _X = aX;
  _Y = aY;
  _Z = aZ;
  _mu = amu;

  R_0 = 1.0;
  //rk = new AdaptSolve();
}
Star::~Star(){
  //delete rk;
}
void Star::Reset(){
  fill(_Dens.begin(), _Dens.end(), 0.0);
  fill(_Temp.begin(), _Temp.end(), 0.0);
  fill(_Mass.begin(), _Mass.end(), 0.0);
  fill(_Lum.begin(), _Lum.end(), 0.0);
  fill(_OptD.begin(), _OptD.end(), 0.0);
  fill(_Rad.begin(), _Rad.end(), 0.0);
  //fill(_Pres.begin(), _Pres.end(), 0.0);
    
  _Kes.clear();
  _KH.clear();
  _Kff.clear();
  _DegPres.clear();
  _GasPres.clear();
  _RadPres.clear();
  _PP.clear();
  _CNO.clear();
  _3a.clear();
  _Pres.clear();
  //_LumDer.clear;
  //_Pder.clear();
  //rk->Reset();
}
void Star::NewStar(double dens, double temp, double aX, double aY, double aZ, double amu){
  central_dens = dens;
  central_temp = temp;
  _X = aX;
  _Y = aY;
  _Z = aZ;
  _mu = amu;
  Reset();
}

double Star::dPdp(double aDens, double aT){//partial der of P wrt density
  double dP = (pow(3.*pow(M_PI,2.),2./3.)/3.)*(pow(hbar,2.))/(me*mp)*pow(aDens/mp,2./3.) + k_b*aT/(_mu*mp);
  return dP;
}

double Star::dPdT(double R,double dens,double temp){
  double dP = dens*k_b/(_mu*mp) + (4./3.)*a*pow(temp,3.0);
  return dP;
  
}

 
//Energy Generation Rates 
//revised
double Star::EGR_PP(double R, double dens, double temp){ // will this function take the density and temperature as vectors or doubles?
  double dens_5 = dens*1e-5;
  double T_6 = temp*1e-6;
  double eps = 1.07e-7*dens_5*pow(_X,2.)*pow(T_6,4.); // not sure how to call X here
  if(isnan(eps) || eps < 1.0e-70){eps = 0.0;}
  _PP.push_back(eps);
  return eps;
}
//revised
double Star::EGR_CNO(double R, double dens, double temp){// same as above fn but for CNO
  double dens_5 = dens*1e-5;
  double T_6 = temp*1e-6;
  double X_cno = 0.03*_X;
  double eps = 8.24e-26*dens_5*_X*X_cno*pow(T_6,19.9);
  if(isnan(eps) || eps < 1.0e-70){eps = 0.0;}
  _CNO.push_back(eps);
  return eps;
}
double Star::EGR_3a(double R, double dens, double temp){// same as above fn but for CNO
  double dens_5 = dens*1e-5;
  double T_8 = temp*1e-8;
  double eps = 3.85e-8*pow(dens_5,2.)*pow(_Y,3.)*pow(T_8,44.0);
  if(isnan(eps) || eps < 1.0e-70){eps = 0.0;}
  _3a.push_back(eps);
  return eps;
  }



//Opacity Function
//revised
double Star::Opacity(double dens, double temp){
  double dens_3 = dens*1e-3;

  double Kes = 0.02*(1+_X);
  double Kff = 1.0e24*(_Z+0.0001)*pow(dens_3,0.7)*pow(temp,-3.5);
  double KH = 2.5e-32*(_Z/0.02)*pow(dens_3,0.5)*pow(temp,9.);
  /*cout << dens << endl;
  cout << ", " << temp;
  cout << ", " << Kes;
  cout << ", " << Kff;
  cout << ", " << pow(dens_3,0.7) << ", " << pow(temp,-3.5);
  cout << ", " << KH << endl;*/
  //_Kes.push_back(Kes);
  //_Kff.push_back(Kff);
  //_KH.push_back(KH);
  double OPsum = pow(KH,-1.) + pow(max(Kes,Kff),-1.);
  return pow(OPsum,-1);
}

void Star::FillOpacity(){
  for(int i = 0; i<_Rad.size();i++){
    double dens_3 = _Dens.at(i)*1.0e-3;;

    double Kes = 0.02*(1+_X);
    double Kff = 1.0e24*(_Z+0.0001)*pow(dens_3,0.7)*pow(_Temp.at(i),-3.5);
    double KH = 2.5e-32*(_Z/0.02)*pow(dens_3,0.5)*pow(_Temp.at(i),9.);
    /*cout << dens << endl;
      cout << ", " << temp;
      cout << ", " << Kes;
      cout << ", " << Kff;
      cout << ", " << pow(dens_3,0.7) << ", " << pow(temp,-3.5);
      cout << ", " << KH << endl;*/
    _Kes.push_back(Kes);
    _Kff.push_back(Kff);
    _KH.push_back(KH);
  }
}

double Star::OpBC(double dens,double temp, double dp){
  double t  = Opacity(dens, temp);
  return t*dens*dens/abs(dp);
}
  
    //Pressure
//revised
double Star::Pressure(double R,double dens,double temp){
  
  //  double P = (pow(3*pow(M_PI,2.),2./3.)/5.)*(pow(hbar,2.)/(me)*pow(dens/mp,5./3.) +
  //					     dens*k_b*temp/(_mu*mp) + (1./3.)*a*pow(temp,4.));
    
  double degp = (pow(3.*M_PI*M_PI, 2./3.)/5.) * (pow(hbar,2.)/(me))*pow(dens/mp,5./3.);
  double gasp = dens*k_b*temp/(_mu*mp);
  double radp = (1./3.)*a*pow(temp,4.);
  return degp+gasp+radp;
}

void Star::FillPres(){
  for(int i = 0; i<_Rad.size();i++){
    double degp = (pow(3.*M_PI*M_PI, 2./3.)/5.) * (pow(hbar,2.)/(me))*pow(_Dens.at(i)/mp,5./3.);
    double gasp = _Dens.at(i)*k_b*_Temp.at(i)/(_mu*mp);
    double radp = (1./3.)*a*pow(_Temp.at(i),4.);
    _DegPres.push_back(degp);
    _GasPres.push_back(gasp);
    _RadPres.push_back(radp);
    _Pres.push_back(degp+gasp+radp);
  }
}

//Derivatives wrt to r
// Mass change with radius
// This also needs to be the density at radius R
double Star::dMdr(double R,double dens){
  double dM = 4.*M_PI*pow(R,2.)*dens;
  //cout << "Mass: " << dM << endl;
  return dM;
}
//revised
double Star::dLdr(double R, double dens, double temp){// luminosity change with radius
  //cout << EGR_CNO(R,dens,temp) << ", " <<  EGR_PP(R,dens,temp) << ", " <<  EGR_3a(R,dens,temp) << endl;
  double dL = 4.*M_PI*pow(R,2.)*dens*(EGR_CNO(R,dens,temp) + EGR_PP(R,dens,temp) + EGR_3a(R,dens,temp));
  //cout << "Lumi: " << dL << endl;
  return dL;
}
//revised
double Star::dtaudr(double dens, double temp){
  double dtau = Opacity(dens,temp)*dens;
  //cout << "OptD: " << dtau << endl;
  return dtau;
}
//revised
double Star::dTdr(double R, double dens, double temp, double mass, double lum){
  /* cout << R;
  cout << ", " << dens;
  cout << ", " << temp;
  cout << ", " << mass;
  cout << ", " << lum;
  cout << ", " << Opacity(dens,temp) << endl;*/
  double rad  = 3.0 *Opacity(dens,temp)*dens*lum / (16*M_PI*a*c*pow(temp,3.)*pow(R,2.));
  double conv = (1. - 1.0/agamma)* temp*G*mass*dens/(Pressure(R,dens,temp)*pow(R,2.));
  //cout << "Rad: " << rad << endl;
  //cout << "Conv: " << conv << endl;
  if(isnan(rad)){throw out_of_range("NaN Temp Grad");}
  return -1.*min(rad,conv);
}
    
// density change with radius
double Star::dpdr(double R,double dens, double temp, double mass,double dt){
  double dp = -1.0*(G*mass*dens*pow(R,-2.) + dPdT(R,dens,temp)*dt)/dPdp(dens,temp);
  //cout << "Dens: " << dp << endl;
  return dp;
}

int Star::MaxArg(){
  return distance(_Rad.begin(),std::max_element(_Rad.begin(),_Rad.end()));
}
// MUST HAVE EVALUATED STAR TO CALL FOLLOWING FUNCTIONS
int Star::SurfRad(){  
  vector<double> dt = _OptD;
  int m = MaxArg();

  //cout << _OptD.size() << endl;
  for(int i = 0; i<m;i++){
    dt.at(i) += abs(( _OptD.at(m)- dt.at(i) - (2./3.)));
    //cout << "Test 1 " << _OptD.at(i) << ", " << dt.at(i) << endl;
  } 			  
  int a = distance(dt.begin(),std::min_element(dt.begin(),dt.begin()+m));
  if(abs(_OptD.at(a)) < 1.e-45){
    a = m;
  }
  return a;
}

double Star::LumBisec(){ 
  int a = SurfRad();
  //cout << "Test 2" << endl;
  double top = _Lum.at(a) - 4.0*M_PI * sigma_sb * pow(_Rad.at(a),2.0)*pow(_Temp.at(a),4.);
  //cout << "Test 3: " << a << ", " << _Rad.at(a) << ", " <<  _Lum.at(a) << ", " << 4.0*M_PI * sigma_sb * pow(_Rad.at(a),2.0)*pow(_Temp.at(a),4.) <<  endl;
  double bot = sqrt(4.0*M_PI*sigma_sb* pow(_Rad.at(a),2.0)*pow(_Temp.at(a),4.)*_Lum.at(a));
  //cout << "Test 3" << endl;
  return top/bot;
}


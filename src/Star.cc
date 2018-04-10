#include "Star.hh"

using namespace std;

Star::Star(double dens, double temp, double aX, double aY, double aZ, double amu){
  central_dens = dens;
  central_temp = temp;
  _X = aX;
  _Y = aY;
  _Z = aZ;
  _mu = amu;

  R_0 = 0.01;
}
Star::Star(const Star &a){
 _Dens = a._Dens;
 _Temp = a._Temp;
 _Mass = a._Mass;
 _Lum = a._Lum;
 _OptD = a._OptD;
 _Rad = a._Rad;
 _dLdr = a._dLdr;
 _dPdT = a._dPdT;
 _dTRad = a._dTRad;
 _dTConv = a._dTConv;

 
 _Kes = a._Kes;
 _KH =  a._KH;
 _Kff = a._Kff;
 _DegPres = a._DegPres;
 _GasPres = a._GasPres;
 _RadPres = a._RadPres;
 _PP = a._PP;
 _CNO = a._CNO;
 _3a = a._3a;
 _Pres = a._Pres;
 _X = a._X;
 _Y = a._Y;
 _Z = a._Z;
 _mu = a._mu;
 central_dens = a.central_dens; 
 central_temp = a.central_temp;
 R_surf = a.R_surf;
 R_0 = a.R_0;
 _MaxRad = a._MaxRad;
}
Star::~Star(){
  //delete rk;
}
void Star::Reset(){
  _Dens.clear();
  _Temp.clear();
  _Mass.clear();
  _Lum.clear();
  _OptD.clear();
  _Rad.clear();
  _dTRad.clear();
  _dTConv.clear();
  _Kes.clear();
  _KH.clear();
  _Kff.clear();
  _DegPres.clear();
  _GasPres.clear();
  _RadPres.clear();
  _Pres.clear();
  _PP.clear();
  _CNO.clear();
  _3a.clear();
  _Pres.clear();
  _dLdr.clear();
  
  R_0 = 1.0;
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
  double dP = (1.0/3.0)*(pow(3.*pow(M_PI,2.),(2.0/3.0))*((pow(hbar,2.0))/(me*mp))*pow((aDens/mp),(2./3.))) + (k_b*aT/(_mu*mp));
  return dP;
}

double Star::dPdT(double R,double dens,double temp){
  double dP = dens*k_b/(_mu*mp) + (4./3.)*a*pow(temp,3.0);
  _dPdT.push_back(dP);
  return dP;
  
}

void Star::FilldPdT(){
    for(int i = 0; i<_Rad.size();i++){
      _dPdT.push_back(_Dens.at(i)*k_b/(_mu*mp) + (4./3.)*a*pow(_Temp.at(i),3.0));
    }
}

void Star::FillEGR(){
    for(int i = 0; i<_Rad.size();i++){
      _PP.push_back(EGR_PP(_Rad.at(i),_Dens.at(i),_Temp.at(i)));
      _CNO.push_back(EGR_CNO(_Rad.at(i),_Dens.at(i),_Temp.at(i)));
      _3a.push_back(EGR_3a(_Rad.at(i),_Dens.at(i),_Temp.at(i)));
    }
}
 
//Energy Generation Rates 
//revised
double Star::EGR_PP(double R, double dens, double temp){ 
  double dens_5 = dens*1.0e-5;
  double T_6 = temp*1.0e-6;
  double eps = 1.07e-7*dens_5*pow(_X,2.)*pow(T_6,4.); 
  if(isnan(eps) || eps < 1.0e-70){eps = 0.0;}
  return eps;
}
//revised
double Star::EGR_CNO(double R, double dens, double temp){// same as above fn but for CNO
  double dens_5 = dens*1.0e-5;
  double T_6 = temp*1.0e-6;
  double X_cno = 0.03*_X;
  double eps = 8.24e-26*dens_5*_X*X_cno*pow(T_6,19.9);
  if(isnan(eps) || eps < 1.0e-70){eps = 0.0;}
  return eps;
}
double Star::EGR_3a(double R, double dens, double temp){// same as above fn but for CNO
  double dens_5 = dens*1.0e-5;
  double T_8 = temp*1.0e-8;
  double eps = 3.85e-8*pow(dens_5,2.)*pow(_Y,3.)*pow(T_8,44.0);
  if(isnan(eps) || eps < 1.0e-70){eps = 0.0;}
  return eps;
  }



//Opacity Function
//revised
double Star::Opacity(double dens, double temp){
  double dens_3 = dens*1.e-3;

  double Kes = 0.02*(1+_X);
  double Kff = 1.0e24*(_Z+0.0001)*pow(dens_3,0.7)*pow(temp,-3.5);
  double KH = 2.5e-32*(_Z/0.02)*pow(dens_3,0.5)*pow(temp,9.);
  if(isnan(Kes) || Kes < 1.0e-15){
    Kes = 1.0e-30;
  }
  if(isnan(Kff) || Kff < 1.0e-15){
    Kff = 1.0e-30;
  }
  if(isnan(KH) || KH < 1.0e-15){
    KH = 1.0e-30;
  }
  double OPsum = pow(KH,-1.) + pow(max(Kes,Kff),-1.);
  return pow(OPsum,-1);
}

void Star::FillOpacity(){
  for(int i = 0; i<_Rad.size();i++){
    double dens_3 = _Dens.at(i)*1.0e-3;;

    double Kes = 0.02*(1+_X);
    double Kff = 1.0e24*(_Z+0.0001)*pow(dens_3,0.7)*pow(_Temp.at(i),-3.5);
    double KH = 2.5e-32*(_Z/0.02)*pow(dens_3,0.5)*pow(_Temp.at(i),9.);
    _Kes.push_back(Kes);
    _Kff.push_back(Kff);
    _KH.push_back(KH);
  }
}

double Star::OpBC(double dens,double temp, double dp){
  double t  = Opacity(dens, temp);
  return t*dens*dens/fabs(dp);
}

//Pressure
//revised
double Star::Pressure(double R,double dens,double temp){
  double degp = (pow(3.*M_PI*M_PI, 2./3.)/5.) * (pow(hbar,2.)/(me))*pow(dens/mp,5./3.);
  double gasp = dens*k_b*temp/(_mu*mp);
  double radp = (1./3.)*a*pow(temp,4.);
  if(isnan(degp) || degp < 1.0e-15){
      degp = 1.0e-30;
    }
  if(isnan(gasp) || gasp < 1.0e-15){
    gasp = 1.0e-30;
  }
  if(isnan(radp) || radp < 1.0e-15){
    radp = 1.0e-30;
  }
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
double Star::dMdr(double R,double dens){
  double dM = 4.*M_PI*pow(R,2.)*dens;
  return dM;
}
//revised
double Star::dLdr(double R, double dens, double temp){// luminosity change with radius
  double dL = 4.*M_PI*pow(R,2.)*dens*(EGR_CNO(R,dens,temp) + EGR_PP(R,dens,temp) + EGR_3a(R,dens,temp));
  return dL;
}

void Star::FilldLdr(){
    for(int i = 0; i<_Rad.size();i++){
      _dLdr.push_back(4.*M_PI*pow(_Rad.at(i),2.)*_Dens.at(i)*(EGR_CNO(_Rad.at(i),_Dens.at(i),_Temp.at(i)) + EGR_PP(_Rad.at(i),_Dens.at(i),_Temp.at(i)) + EGR_3a(_Rad.at(i),_Dens.at(i),_Temp.at(i))));
    }
}
//revised
double Star::dtaudr(double dens, double temp){
  double dtau = Opacity(dens,temp)*dens;
  if(dens > 100.0){
    return 0 ;
  } else{
    return dtau;
  }
}
//revised
double Star::dTdr(double R, double dens, double temp, double mass, double lum){
  double rad  = 3.0 *Opacity(dens,temp)*dens*lum / (64.*M_PI*sigma_sb*pow(temp,3.)*pow(R,2.));
  double conv = (1. - pow(agamma,-1.0))*temp*G*mass*dens/(Pressure(R,dens,temp)*pow(R,2.));
  if(isnan(rad)){throw out_of_range("NaN Temp Grad");}
  return -1.0*min(rad,conv);
}

void Star::FilldTdr(){
   for(int i = 0; i<_Rad.size();i++){
     _dTRad.push_back(3.0 *Opacity(_Dens.at(i),_Temp.at(i))*_Dens.at(i)*_Lum.at(i) / (64.*M_PI*sigma_sb*pow(_Temp.at(i),3.)*pow(_Rad.at(i),2.)));
     _dTConv.push_back((1. - pow(agamma,-1.0))*_Temp.at(i)*G*_Mass.at(i)*_Dens.at(i)/(Pressure(_Rad.at(i),_Dens.at(i),_Temp.at(i))*pow(_Rad.at(i),2.)));
   }
}
    
// density change with radius
double Star::dpdr(double R,double dens, double temp, double mass,double dt){
  double dp = -1.0*(G*mass*dens/pow(R,2.) + dPdT(R,dens,temp)*dt)/dPdp(dens,temp);
  return dp;
}




int Star::MaxArg(){
  return distance(_Rad.begin(),std::max_element(_Rad.begin(),_Rad.end()));
}
// MUST HAVE EVALUATED STAR TO CALL FOLLOWING FUNCTIONS
int Star::SurfRad(){    
  int m = MaxArg();
  vector<double> dt;
  int a =0;
  for(int i = 0; i<=m;i++){
    dt.push_back(abs((_OptD.at(m)-_OptD.at(i)  - (2.0/3.0))));
  }
  
  double num = *std::min_element(dt.begin(),dt.end());
  a = distance(dt.begin(),std::min_element(dt.begin(),dt.end()));
  dt.clear();
  _MaxRad = a;
  return a;
}

double Star::LumBisec(){ 
  int a = SurfRad();
  double top = _Lum.at(a) - (4.0*M_PI*sigma_sb*pow(_Rad.at(a),2.0)*pow(_Temp.at(a),4.));
  double bot = sqrt(4.0*M_PI*sigma_sb* pow(_Rad.at(a),2.0)*pow(_Temp.at(a),4.)*_Lum.at(a));
  return top/bot;
}


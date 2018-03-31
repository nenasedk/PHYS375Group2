/*
 * Constants Class
 *
 * This class will define all the units used in the Stars project,
 * and will allow for easy unit conversion.
 *
 */

#ifndef _Constants_
#define _Constants_

using namespace std;

// If we want to get a bit fancier. Probably unnecssary.
//struct Unit{long double value; string unit;};
  
// Units, Constants, and conversions in C++
// Base Sun Units
static const long double M_Sun = 1.988435e30; // kg
//static const long double M_Sun = 1.0; // MSun
static const long double Ls = 3.848e26; // W
//static const long double Ls = 1.0; // LSun
static const long double Tes = 5778.;    // K

//elementary particles
static const long double e = 1.602e-19; // C elementary charge
static const long double mp = 1.673e-27;// kg
static const long double me = 9.109e-31;// kg
  
//Distances (in km)
static const long double AU = 1.496e8; // m
static const long double pc = 3.086e13;// m

//Other Constants
static const long double  G = 6.673e-11; //
static const long double  eps_0 = 8.854e-12; //permittivity of vacuum
static const long double  mun_0 =(3.14159)*4e-7;//permeability of vacuum
static const long double  c = 2.998e8; // m/s
static const long double  h = 6.626e-34; // Js
static const long double  hbar = 1.055e-34; // Js
static const long double  k_b = 1.381e-23; // J/k
static const long double  sigma_sb = 5.670e-8; //stefan Boltzman constant
const long double  a = 4.*sigma_sb/c; // radiative constant

//Things from Assignment
//XCNO = 0.03*X;
//ro = 1 #temporary placeholder
//ro5 = ro*1e-5
//ro3 = ro*1e-3
//T = 1 #temporary placeholder
//T6 = T*1e-6
//ePP = (1.07e-7)*(ro5)*(X**2)*(T6**4)
//eCNO = (8.24e-26)*(ro5)*X*XCNO*(T6**19.9)
//eroT = ePP + eCNO
//a = 4*ksb/c #radiative pressure coefficient
static const long double agamma = 5./3.; // adiabatic constant

#endif

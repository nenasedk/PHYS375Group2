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
//struct Unit{double value; string unit;};
  
// Units, Constants, and conversions in C++
// Base Sun Units
static const double M_Sun = 1.988435e30; // kg
static const double Ls = 3.848e26; // W
static const double Tes = 5778;    // K

//elementary particles
static const double e = 1.602e-19; // C elementary charge
static const double mp = 1.673e-27;// kg
static const double me = 9.109e-31;// kg
  
//Distances (in km)
static const double AU = 1.496e8; // m
static const double pc = 3.086e13;// m

//Other Constants
static const double  G = 6.673e-11; //
static const double  eps = 8.854e-12; //permittivity of vacuum
static const double  mun =(3.14159)*4e-7;//permeability of vacuum
static const double  c = 2.998e8; // m/s
static const double  h = 6.626e-34; // Js
static const double  hbar = 1.055e-34; // Js
static const double  k = 1.381e-23; // J/k
static const double  sigma_sb = 5.670e-8; //stefan Boltzman constant
const double  a = 4*sigma_sb/c; // radiative constant

//Mean Molecular Weight (set for sun)
//X = 0.734;
//Y = 0.250;
//Z = 0.016;
//mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1);

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
static const double agamma = 5/3; // adiabatic constant

//Rosseland Mean Opacities
//kapes = 0.02*(1+X)
//kapff = (1.0e24)*(Z + 0.0001)*(ro3**0.7)*(T**(-3.5))
//kapH = (2.5e-32)*(Z/0.02)*(ro3**0.5)*T**9

#endif

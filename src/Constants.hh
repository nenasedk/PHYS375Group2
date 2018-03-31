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
//static const double M_Sun = 1.0; // MSun
static const double Ls = 3.848e26; // W
//static const double Ls = 1.0; // LSun
static const double Tes = 5778.;    // K

//elementary particles
static const double e = 1.602e-19; // C elementary charge
static const double mp = 1.673e-27;// kg
static const double me = 9.109e-31;// kg
  
//Distances (in km)
static const double AU = 1.496e8; // m
static const double pc = 3.086e13;// m

//Other Constants
static const double  G = 6.673e-11; //
static const double  eps_0 = 8.854e-12; //permittivity of vacuum
static const double  mun_0 =(3.14159)*4e-7;//permeability of vacuum
static const double  c = 2.998e8; // m/s
static const double  h = 6.626e-34; // Js
static const double  hbar = 1.055e-34; // Js
static const double  k_b = 1.381e-23; // J/k
static const double  sigma_sb = 5.670e-8; //stefan Boltzman constant
const double  a = 4.*sigma_sb/c; // radiative constant
static const double agamma = 5./3.; // adiabatic constant

#endif

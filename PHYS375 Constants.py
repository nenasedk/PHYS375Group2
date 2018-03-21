##Units, Constants, and conversions in Python
from astropy import units as u
from astropy import constants as const

# The Astropy constants library contains many useful
# constants, complete with units.
#
# The units library allows quick conversion between
# quantities with the same dimensions.

#Base Sun Units

# Example for using astropy
Ms = const.M_sun       # By using const we can use the astropy values 
Mskg = Ms.to(u.M_sun)  # units allows us to convert from kg to solar masses
Ls = const.L_sun
Tes = 5778*u.K         # You can multiply units onto any numerical quantity

#elementary particles
e = 1.602e-19 #elementary charge
mp = 1.673e-27
me = 9.109e-31

#Distances (in km)
AU = 1.496e8
pc = 3.086e13

#Other Constants
G = 6.673e-11
eps = 8.854e-12 #permittivity of vacuum
mun =(3.14159)*4e-7 #permeability of vacuum
c = 2.998e8
h = 6.626e-34
hbar = 1.055e-34
k = 1.381e-23
ksb = 5.670e-8 #stefan Boltzman constant

#Mean Molecular Weight (set for sun)
X = 0.734
Y = 0.250
Z = 0.016
mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1)

#Things from Assignment
XCNO = 0.03*X
ro = 1 #temporary placeholder
ro5 = ro*1e-5
ro3 = ro*1e-3
T = 1 #temporary placeholder
T6 = T*1e-6
ePP = (1.07e-7)*(ro5)*(X**2)*(T6**4)
eCNO = (8.24e-26)*(ro5)*X*XCNO*(T6**19.9)
eroT = ePP + eCNO
a = 4*ksb/c #radiative pressure coefficient

#Rosseland Mean Opacities
kapes = 0.02*(1+X)
kapff = (1.0e24)*(Z + 0.0001)*(ro3**0.7)*(T**(-3.5))
kapH = (2.5e-32)*(Z/0.02)*(ro3**0.5)*T**9





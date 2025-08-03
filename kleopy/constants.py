"""
Constants (:mod:`kleopy.constants`)
==========================================================
Physical constants and parameters used in the kleopy package.

"""
from numpy import pi

#----- Universal Constants -----
G = 6.67430e-11  #Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    #Speed of light in vacuum (m/s)

#----- 216-Kleopatra asteroid parameters -----
m1 = 1.1014e+18  #Mass of first body in 216 Kleopatra (kg)
m2 = 1.0350e+18  #Mass of second body in 216 Kleopatra (kg)
ms = 4.1547e+17  #Mass of 216 Kleopatra segment (kg)
l = 1            #Length of 216 Kleopatra segment (dimensionless, normalized to l)
l1 = 0.486608    #Distance of first body from COM in 216 Kleopatra (dimensionless, normalized to l)
l2 = 0.513392    #Distance of second body from COM in 216 Kleopatra (dimensionless, normalized to l)
T = 2 * pi       #Rotation period of 216 Kleopatra (dimensionless, normalized to 1/omega or T/2pi)
kappa = 0.991    #Kappa parameter for 216 Kleopatra (dimensionless) ##Check again later
mu = 0.484       #Measure of how imbalanced the mass distribution is (towards 2nd body) (dimensionless)
mu_s = 0.163     #Measure of how much the segment contributes to the total mass (dimensionless)

#----- Physical values -----
#In the paper, the lengths and time are nondimensional and uses l and 1/omega as units, respectively.
#Thus, for converting convenience, we provide the literature values in meters and seconds.
l1_lit = 57322.422399999996     #Distance of first body from COM in 216 Kleopatra (m)
l2_lit = 60477.5776             #Distance of second body from COM in 216 Kleopatra (m)
l_lit = 117800                  #Length of 216 Kleopatra segment (m)
T_lit = 5.385                   #Rotation period of 216 Kleopatra (hours)
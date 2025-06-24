"""
Constants (:mod:`kleopy.constants`)
==========================================================
Physical constants and parameters used in the kleopy package.

"""

#----- Universal Constants -----
G = 6.67430e-11  #Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    #Speed of light in vacuum (m/s)

#----- 216-Kleopatra asteroid parameters -----
m1 = 1.1014e+18  #Mass of first body in 216 Kleopatra (kg)
m2 = 1.0350e+18  #Mass of second body in 216 Kleopatra (kg)
ms = 4.1547e+17  #Mass of 216 Kleopatra segment (kg)
l = 117800       #Length of 216 Kleopatra segment (m)
l1 = 57322.422399999996    #Distance of first body from COM in 216 Kleopatra (m)
l2 = 60477.5776  #Distance of second body from COM in 216 Kleopatra (m)
kappa = 0.991    #Kappa parameter for 216 Kleopatra (dimensionless) ##Check again later
T = 5.385        #Rotation period of 216 Kleopatra (hours)
mu = 0.484       #Measure of how imbalanced the mass distribution is (towards 2nd body) (dimensionless)
mu_s = 0.163     #Measure of how much the segment contributes to the total mass (dimensionless)

#----- Literature values -----
#I noticed that l1 and l2 adds up to 1, that may indicate that they are fractions of the total length l
#thus i will put them here for reference
l1_lit = 0.486608
l2_lit = 0.513392
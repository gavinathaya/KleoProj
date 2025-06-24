"""
Orbital equations (:mod:`kleopy.orbital_eq`)
==========================================================
Various equations related to the orbital dynamics of 216-Kleopatra.

Includes potentials and equations of motion.
"""
from kleopy.constants import (
    G, m1, m2, ms, l, l1, l2, kappa, T, mu, mu_s
)
import numpy as np

#Bound all constants to the module to improve performance
G = G; m1 = m1; m2 = m2; ms = ms; l = l; l1 = l1; l2 = l2
kappa = kappa; T = T; mu = mu; mu_s = mu_s

M = m1 + m2 + ms  #Total mass of 216-Kleopatra
from scipy.differentiate import derivative
#----- Potentials -----
def potential(x, y, z) -> float:
    """
    Gravitational potential from 216-Kleopatra at position (x, y, z).

    Parameters
    ----------
    x : float
        x-coordinate in the synodic frame.
    y : float
        y-coordinate in the synodic frame.
    z : float
        z-coordinate in the synodic frame.
    
    Returns
    -------
    U : float
        Gravitational potential at (x,y,z).
    """
    #Calculating distances from the two bodies in 216-Kleopatra
    r1 = np.sqrt((x + l1)**2 + y**2 + z**2) #Distance to first body
    r2 = np.sqrt((x - l2)**2 + y**2 + z**2) #Distance to second body

    #note: i'm 95% sure the "log" in the paper means ln, not base 10...
    # ... because classical gravity is analogous to electrostatics,
    # ... and there the potential of a cylinder is proportional to ln(r),
    # ... so we use np.log which is ln
    U = -G * M * ((1 - mu) * (1-mu_s)/r1 + mu * (1-mu_s)/r2 + mu_s/l * np.log((r1+r2+l)/(r1+r2-l)))
    return U

def potential_eff(x, y, z) -> float:
    """
    Effective gravitational potential from 216-Kleopatra at position (x, y, z) due to synodic frame.

    Parameters
    ----------
    x : float
        x-coordinate in the synodic frame.
    y : float
        y-coordinate in the synodic frame.
    z : float
        z-coordinate in the synodic frame.
    
    Returns
    -------
    Omega : float
        Effective gravitational potential at (x,y,z).
    """
    #Calculating distances from the two bodies in 216-Kleopatra
    r1 = np.sqrt((x + l1)**2 + y**2 + z**2) #Distance to first body
    r2 = np.sqrt((x - l2)**2 + y**2 + z**2) #Distance to second body

    #note: i'm 95% sure the "log" in the paper means ln, not base 10...
    # ... because classical gravity is analogous to electrostatics,
    # ... and there the potential of a cylinder is proportional to ln(r),
    # ... so we use np.log which is ln
    Omega = 1/2 * (x**2 + y**2) + kappa * ((1 - mu) * (1-mu_s)/r1 + mu * (1-mu_s)/r2 + mu_s/l * np.log((r1+r2+l)/(r1+r2-l)))
    return Omega


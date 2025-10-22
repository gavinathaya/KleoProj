"""
Integration tools (:mod:`kleopy.integrator`)
==========================================================
Various functions related to differential equation integrators for 216-Kleopatra.


Functions
---------

find_dyt(x0, C) -> np.ndarray       : Calculate the dy0t value for a given x0 and C (Jacobian constant).
init_grid(x0_min, x0_max,
        C_min, C_max, 
        dif_x0, dif_C, 
        retgrid=False
        ) -> tuple[np.ndarray, ...] : Initializes a grids of x0, C, and dy0t values for the grid search process.
"""

#Import libraries
import numpy as np; import numpy.typing as npt
from scipy.integrate import solve_ivp

def find_dxt(fun, x0: npt.ArrayLike, dy0t: npt.ArrayLike, method: str = 'DOP853') -> np.ndarray:
    """
    """
    sol = solve_ivp(fun, (0, np.pi), [x0, 0, 0, 0, dy0t, 0], method=method, t_eval = (np.pi,))
    return sol.y[3,0]
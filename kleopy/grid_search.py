"""
Grid Search tools (:mod:`kleopy.grid_search`)
==========================================================
Various equations and functions related to the grid search process of 216-Kleopatra.


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
from scipy.optimize import brentq
from psutil import cpu_count
import sys; import time
from itertools import repeat
from kleopy.orbital_eq import potential_eff, EOM
from kleopy.misc import progress_bar

def find_dy0t(x0: float | npt.ArrayLike, C: float | npt.ArrayLike) -> np.ndarray:
    """
    Calculate the dy0t value for a given x0 and C (Jacobian constant).
    Uses vectorized operations for performance.

    Parameters
    ----------
    x0 : float or np.ndarray
        Initial x-coordinate in the synodic frame.
    C : float or np.ndarray
        Jacobian constant.
    
    Returns
    -------
    dy0t : np.ndarray
        The dy0t value calculated from the effective potential at (x0, 0, 0) and the Jacobian constant C.
    """
    #Ignore invalid values in potential_eff
    np.seterr(invalid='ignore', divide='ignore')

    #Insure inputs are numpy arrays
    x0 = np.asarray(x0)
    C = np.asarray(C)

    #Invalid value handling
    sqrtval = 2 * potential_eff(x0, 0, 0) - C
    badmask = (sqrtval < 0) | ~np.isfinite(sqrtval) #Mask for invalid values
    
    #Calculate dy0t
    dy0t = np.sqrt(sqrtval)

    #Set invalid values to NaN
    dy0t[badmask] = np.nan 
    return dy0t

def find_C(x0: float | npt.ArrayLike, dy0t: float | npt.ArrayLike) -> np.ndarray:
    """
    Calculate the dy0t value for a given x0 and C (Jacobian constant).
    Uses vectorized operations for performance.

    Parameters
    ----------
    x0 : float or np.ndarray
        Initial x-coordinate in the synodic frame.
    dy0t : float or np.ndarray
        Initial dyt constant.
    
    Returns
    -------
    C : np.ndarray
        The Jacobian constant calculated from the effective potential at (x0, 0, 0) and the initial dyt.
    """
    #Ignore invalid values in potential_eff
    np.seterr(invalid='ignore', divide='ignore')

    #Insure inputs are numpy arrays
    x0 = np.asarray(x0)
    dy0t = np.asarray(dy0t)

    #Invalid value handling
    C = 2 * potential_eff(x0, 0, 0) - dy0t**2
    badmask = ~np.isfinite(C) #Mask for invalid values

    #Set invalid values to NaN
    C[badmask] = np.nan 
    return C

def old_init_grid(x0_min: float = -3, x0_max: float = 2, C_min: float = -3, C_max: float = 5, *, 
              dif_x0: float = 0.001, dif_C:float = 0.001, 
              retgrid = False) -> tuple[np.ndarray, ...]:
    """
    Initializes a grids of x0, C, and dy0t values for the grid search process.
    Indexing is done in ij format, meaning the i-th x0 and j-th C.
    Uses vectorized operations for performance.

    Parameters
    ----------
    x0_min : float, optional
        Minimum value for x0. Default is -3.
    x0_max : float, optional
        Maximum value for x0. Default is 2.
    C_min : float, optional
        Minimum value for C. Default is -3.
    C_max : float, optional
        Maximum value for C. Default is 5.
    dif_x0 : float, optional
        Step size for x0. Default is 0.001.
    dif_C : float, optional
        Step size for C. Default is 0.001.

    Returns
    -------
    X0_grid : np.ndarray
        2D array of shape (nx0, nC) containing x0 values for i-th x0 and j-th C.
        Where C is the jacobian constant.
    C_grid : np.ndarray
        2D array of shape (nx0, nC) containing C values for i-th x0 and j-th C.
        Where C is the jacobian constant.
    DY0T_grid : np.ndarray
        2D array of shape (nx0, nC) containing dy0t values for i-th x0 and j-th C.
        Where C is the jacobian constant.
    grid : np.ndarray, optional
        Only returns if retgrid is true.
        Array of shape (N, 2) containing pairs of (x0, C).
    """
    #Start time and confirm grid initialization start
    start_time = time.time()
    sys.stdout.write(f"{'\033[94m'}Initializing grid...{'\033[0m'}")
    #Create a grid of x0 and C values as inputs for the find_dyt function
    x0_array = np.linspace(x0_min, x0_max, int((x0_max - x0_min) / dif_x0) + 1)
    C_array = np.linspace(C_min, C_max, int((C_max - C_min) / dif_C) + 1)
    X0_grid, C_grid = np.meshgrid(x0_array, C_array, indexing='ij') #ij indexing: i-th x0, j-th C
    
    #Apply the find_dyt function to each pair of (x0, C) in the grid
    DY0T_grid = find_dy0t(X0_grid, C_grid)

    #End time and confirm grid initialization
    end_time = time.time()
    sys.stdout.write(f"\r{'\033[92m'}Grid initialized! Completed in {end_time - start_time} seconds{'\033[0m'} ")

    #Return the grids and optionally the grid of pairs (x0, C)
    if retgrid is True:
        # Array of shape (N, 2) containing pairs of (x0, C).
        grid = np.array(np.meshgrid(x0_array, C_array)).T.reshape(-1, 2)
        return X0_grid, C_grid, DY0T_grid, grid
    else:
        #Returns ij grid of x0, C, and dy0t values.
        return X0_grid, C_grid, DY0T_grid

def _old_process_grid(fun, X0_grid: np.ndarray, DY0T_grid: np.ndarray): #-> np.ndarray:
    """
    Integrates the grid of x0, C, and dy0t values until the T/2 point.
    Uses vectorized operations for performance.

    Parameters
    ----------
    X0_grid : np.ndarray
        2D array of shape (nx0, nC) containing x0 values for i-th x0 and j-th C.
        Where C is the jacobian constant.

    DY0T_grid : np.ndarray
        2D array of shape (nx0, nC) containing dy0t values for i-th x0 and j-th C.
        Where C is the jacobian constant.
    
    Returns
    -------
    DXT_grid : np.ndarray
        2D array of shape (nx0, nC) containing dxt(T/2) values for i-th x0 and j-th C.
        Where C is the jacobian constant.
    """
    #Insure inputs are numpy arrays
    X0_grid = np.asarray(X0_grid)
    DY0T_grid = np.asarray(DY0T_grid)

    #Calculate dxt(T/2) for each pair of (x0, C) in the grid
    DXT_grid = np.zeros_like(DY0T_grid)  # Initialize the grid
    

    return DXT_grid

def EOM_event(t,Y):
        return Y[1]
EOM_event.terminal = True # pyright: ignore[reportFunctionMemberAccess]
EOM_event.direction = -1  # pyright: ignore[reportFunctionMemberAccess]

def _grid_search_worker(x0_array: np.ndarray, dy0t_array: np.ndarray, *, time_span = (0, 1e6), event_func = EOM_event, atol = 1e-12, rtol = 1e-10, debug = False):
    """
    Single worker to do the grid search process
    """
    # TODO: Dimension check, both should be 1 dimensional NDArrays and have the same size
    i_total = x0_array.size
    found_orbits = []
    #Python list for storing dxt at T/2 values
    dxt_at_half_T = []
    start_time = time.time()
    for i, x0_val in enumerate(x0_array):
        dy0t = dy0t_array[i]
        #Check if dy0t is Nan
        if dy0t == np.nan:
            dxt_at_half_T.append(np.nan)
            continue #Continue to the next x0 value
        Y0 = [x0_val, 0, 0, dy0t]
        sol = solve_ivp(EOM, time_span, Y0, method = 'Radau', events = event_func, atol = atol, rtol = rtol)
        if sol.y_events[0].size < 1:
            print(f"{'\033[31m'}WARNING: y = 0 NEVER REACHED (REALLY STIFF) FOR INITIAL CONDITIONS x0: {x0_val}| dy0t: {dy0t}{'\033[0m'}")
            dxt_at_half_T.append(np.nan)
            continue
        dxt_found = sol.y_events[0][1] #0 index as the "first" event function triggered and 1 as the "second" component of Y state vector
        dxt_at_half_T.append(dxt_found)
        current_time = time.time()
        progress_bar(i,i_total)
        sys.stdout.write(f"Elapsed time: {current_time - start_time}")
    dxt_at_half_T = np.asarray(dxt_at_half_T)
    sign_changes = np.where(dxt_at_half_T[:-1] * dxt_at_half_T[1:] < 0)[0] # pyright: ignore[reportIndexIssue]
    exact_zeros = np.where(dxt_at_half_T == 0)[0]
    if exact_zeros.size > 0:
        found_orbits.extend(zip(x0_array[exact_zeros], dy0t_array[exact_zeros]))
    
    return None

def process_grid(x0_min: float = -3, x0_max: float = 2, C_min: float = -3, C_max: float = 5, *, 
              dif_x0: float = 0.001, dif_C:float = 0.001, 
              retgrid = False):
    """
    Initializes a grids of x0, C, and dy0t values for the grid search process.
    Indexing is done in ij format, meaning the i-th x0 and j-th C.
    Uses vectorized operations for performance.

    Parameters
    ----------
    x0_min : float, optional
        Minimum value for x0. Default is -3.
    x0_max : float, optional
        Maximum value for x0. Default is 2.
    C_min : float, optional
        Minimum value for C. Default is -3.
    C_max : float, optional
        Maximum value for C. Default is 5.
    dif_x0 : float, optional
        Step size for x0. Default is 0.001.
    dif_C : float, optional
        Step size for C. Default is 0.001.

    Returns
    -------
    X0_grid : np.ndarray
        2D array of shape (nx0, nC) containing x0 values for i-th x0 and j-th C.
        Where C is the jacobian constant.
    C_grid : np.ndarray
        2D array of shape (nx0, nC) containing C values for i-th x0 and j-th C.
        Where C is the jacobian constant.
    DY0T_grid : np.ndarray
        2D array of shape (nx0, nC) containing dy0t values for i-th x0 and j-th C.
        Where C is the jacobian constant.
    grid : np.ndarray, optional
        Only returns if retgrid is true.
        Array of shape (N, 2) containing pairs of (x0, C).
    """
    #Start time and confirm grid initialization start
    start_time = time.time()
    sys.stdout.write(f"{'\033[94m'}Initializing grid...{'\033[0m'}")
    #Create a grid of x0 and C values
    x0_array = np.linspace(x0_min, x0_max, int((x0_max - x0_min) / dif_x0) + 1)
    C_array = np.linspace(C_min, C_max, int((C_max - C_min) / dif_C) + 1)
    X0_grid, C_grid = np.meshgrid(x0_array, C_array, indexing='ij') #ij indexing: i-th x0, j-th C
    
    #Apply the find_dyt function to each pair of (x0, C) in the grid
    DY0T_grid = find_dy0t(X0_grid, C_grid)

    #End time and confirm grid initialization
    end_time = time.time()
    sys.stdout.write(f"\r{'\033[92m'}Grid initialized! Completed in {end_time - start_time} seconds{'\033[0m'} ")
    return   
import os
import numpy as np
from tqdm import tqdm
import rasterio

from rasterio.mask import mask, raster_geometry_mask
import time

''' 
This file contains all the necessary functions to compute the accessibility metrics
'''

def shortest_travel_time(tt_rast_path):
    """
    Computes the shortest travel time across overlapping travel time rasters.

    This function iterates through all `.tif` or `.geotiff` files in the specified directory, 
    reading travel time raster data and calculating the minimum travel time for each grid cell 
    where the rasters overlap. The result is a raster representing the shortest travel time 
    for each location.

    Parameters:
    ----------
    tt_rast_path : str
        The directory containing `.tif` or `.geotiff` raster files representing travel time data.

    Returns:
    -------
    np.ndarray
        A 2D numpy array where each cell value represents the minimum travel time across all 
        raster files for the corresponding location.
    """
    i = 0
    for file in tqdm(os.listdir(tt_rast_path), desc='Computing shortest travel time'):

        if file.endswith('.tif') or file.endswith('.geotiff'):
            point_rast = rasterio.open(tt_rast_path+file)
            # resample data to target shape
            point_tt_arr = point_rast.read()[0]
            if i == 0:
                min_raster = point_tt_arr.copy()
            else:
                min_raster = np.fmin(point_tt_arr, min_raster)
            i += 1
        else:
            continue
                
        # if i % 100 == 0:
        #     print(f"{i}/{len(os.listdir(tt_rast_path))} Processed rasters")
    return min_raster

#impedance functions for accessibility metrics

def exp_impedance_fun(dist_arr, beta):
    return np.exp(-beta * dist_arr)

def gauss_impedance_fun(dist_arr, beta):
    return np.exp(-beta * dist_arr**2)

# def compute_grav_acc(beta, pop_arr, tt_dir):
#     """
#     Computes gravity-based accessibility weighed by population for each pixel.

#     This function processes travel time raster files in the specified directory and computes 
#     gravity-based accessibility for each grid cell using an exponential impedance function. 
#     The accessibility to each facility is weighed by the total population-weighted .

#     Parameters:
#     ----------
#     beta : float
#         The impedance parameter for the exponential decay function, which controls how accessibility 
#         decreases with increasing travel time.
#     pop_arr : np.ndarray
#         A 2D numpy array representing the population distribution. Cells with non-positive values 
#         or non-finite values are treated as unpopulated and excluded from calculations.
#     tt_dir : str
#         The directory containing `.tif` raster files representing travel time data for different facilities.

#     Returns:
#     -------
#     np.ndarray
#         A 2D numpy array representing the gravity-based accessibility matrix weighed by population, where each cell value 
#         reflects the cumulative normalized accessibility across all facilities.
#     """

    
#     access_mat = np.zeros_like(pop_arr)  # matrix of a_i
    
#     # Precompute the population mask to avoid redundant operations
#     pop_arr_cop = np.where(pop_arr < 0, 0, pop_arr.copy())
#     pop_arr_cop[~np.isfinite(pop_arr_cop)] = 0
    
#     for i, file in enumerate(os.listdir(tt_dir)):
#         if not file.endswith('.tif'):
#             continue
            
#         dist_rast = rasterio.open(tt_dir+file)
#         dist_arr = dist_rast.read()[0]
#         acc_array = exp_impedance_fun(dist_arr, beta)
        
#         # Apply population mask to accessibility
#         acc_array[pop_arr_cop == 0] = 0
#         acc_array[~np.isfinite(acc_array)] = 0
#         denom = np.sum(pop_arr_cop*acc_array)
#         acc_array = acc_array/denom
        
#         # User-based accessibility (cumulative)
#         access_mat += acc_array

#         if i % 100 == 0:
#             print(f"Processed facility {i}/{len(os.listdir(tt_dir))}")
#     access_mat[pop_arr < 0] = np.nan
#     return access_mat
    
# def compute_grav_acc_cumul(beta, pop_arr, tt_dir):
#     """
#     Computes cumulative gravity-based accessibility using an exponential impedance function.

#     This function processes travel time raster files in the specified directory and computes 
#     cumulative accessibility for each grid cell based on an exponential decay function. 
#     The result is a cumulative accessibility matrix where accessibility to all facilities 
#     is aggregated for each grid cell.

#     Parameters:
#     ----------
#     beta : float
#         The impedance parameter for the exponential decay function, which controls how accessibility 
#         decreases with increasing travel time.
#     pop_arr : np.ndarray
#         A 2D numpy array representing the population distribution. Cells with non-positive values 
#         are treated as unpopulated and excluded from calculations.
#     tt_dir : str
#         The directory containing `.tif` raster files representing travel time data for different facilities.

#     Returns:
#     -------
#     np.ndarray
#         A 2D numpy array representing the cumulative gravity-based accessibility matrix, where each cell value 
#         reflects the total accessibility to all facilities based on the exponential decay function.
#     """
    
#     access_mat = np.zeros_like(pop_arr)  # matrix of a_i
    
#     # Precompute the population mask to avoid redundant operations
#     pop_arr_cop = np.where(pop_arr < 0, 0, pop_arr.copy())
    
#     # Precompute the square of the coordinates to avoid repetitive calculation

#     for i, file in enumerate(os.listdir(tt_dir)):
#         if not file.endswith('.tif'):
#             continue

#         dist_rast = rasterio.open(tt_dir+file)
#         dist_arr = dist_rast.read()[0]

#         # Compute accessibility using vectorized operations
#         acc_array = exp_impedance_fun(dist_arr, beta)
        
#         # Apply population mask to accessibility
#         acc_array[pop_arr_cop == 0] = 0
#         acc_array[~np.isfinite(acc_array)] = 0

#         # User-based accessibility (cumulative)
#         access_mat += acc_array

#         if i % 100 == 0:
#             print(f"Processed facility {i}/{len(os.listdir(tt_dir))}")
#     access_mat[pop_arr < 0] = np.nan
#     return access_mat

def compute_cumulative(max_tt:float, mask_arr, tt_dir:str):
    '''
    Computes the cumulative accessibility matrix for markets within a given maximum travel time.

    This function iterates through all `.tif` files in the specified directory and computes 
    how many markets (represented as raster files) are within `max_tt` minutes of travel time 
    for each grid cell. The result is a cumulative accessibility matrix.

    Parameters:
    ----------
    max_tt : float
        The maximum travel time (in minutes) to consider when determining accessibility.
    tt_dir : str
        The directory containing `.tif` raster files representing travel time data for different markets.
    mask_arr : np.ndarray
        Friction surface to use as mask for the country border
    Returns:
    -------
    np.ndarray
        A 2D numpy array representing the cumulative accessibility matrix, where each cell value 
        indicates the number of markets accessible within `max_tt` minutes.
    '''
    i = 0
    for file in tqdm(os.listdir(tt_dir), desc = "Computing cumulative acc."):
        if not file.endswith('.tif') or file.endswith('.geotiff'):
            continue

        dist_rast = rasterio.open(tt_dir+file)
        dist_arr = dist_rast.read()[0]
        if i == 0: 
            access_mat = np.zeros_like(dist_arr)
        # Compute accessibility using vectorized operations
        acc_array = np.zeros_like(dist_arr)
        acc_array[dist_arr<=max_tt] = 1
        acc_array[~np.isfinite(acc_array)] = 0

        # User-based accessibility (cumulative)

        access_mat += acc_array

        # if i % 100 == 0:
        #     print(f"Processed facility {i}/{len(os.listdir(tt_dir))}")
        i += 1
    # set places where the friction surface is infinite to nan
    access_mat[~np.isfinite(mask_arr)] = np.nan
    return access_mat

# def compute_entropy_acc(beta, mask_arr, tt_dir):
#     """
#     Computes entropy-based accessibility (S) and Average Travel-Time (internal energy) (U) for a given population 
#     and travel time rasters.

#     This function processes travel time raster files in the specified directory, calculates impedance-based 
#     accessibility for each facility, and computes two metrics:
#     - `S` (Entropy-based accessibility): Measures diversity in accessibility across multiple facilities.
#     - `U` (Average Travel Time): Measures the average travel time across all facilities according to travel behaviour.

#     Parameters:
#     ----------
#     beta : float
#         The impedance parameter for the exponential decay function, which controls how accessibility 
#         decreases with increasing travel time.
#     mask_arr : np.ndarray
#         Friction surface to use as mask for the country border
#     tt_dir : str
#         The directory containing `.tif` raster files representing travel time data for different facilities.

#     Returns:
#     -------
#     tuple (np.ndarray, np.ndarray)
#         S : np.ndarray
#             A 2D numpy array representing entropy-based accessibility for each grid cell.
#         U : np.ndarray
#             A 2D numpy array representing average travel time for each grid cell.
#     """
    
#     # Precompute the population mask to avoid redundant operations
#     mask_arr_cop = np.where(mask_arr < 0, np.inf, mask_arr.copy())
    
#     # Precompute the square of the coordinates to avoid repetitive calculation
    
#     Z = np.zeros_like(mask_arr_cop)
#     U = np.zeros_like(mask_arr_cop)

#     for i, file in enumerate(os.listdir(tt_dir)):
#         if not file.endswith('.tif'):
#             continue

#         dist_rast = rasterio.open(tt_dir+file)
#         dist_arr = dist_rast.read()[0]

#         # Compute accessibility using vectorized operations
#         acc_array = exp_impedance_fun(dist_arr, beta)
#         acc_array = np.nan_to_num(acc_array, nan = 0)
#         # set places where the friction surface is infinite to nan
#         acc_array[~np.isfinite(mask_arr)] = np.nan
#         #dealing with nans
#         acc_array[~np.isfinite(acc_array)] = 0

#         dist_arr = np.where(dist_arr == np.inf, 0, dist_arr)

#         U += acc_array*dist_arr
        
#         Z = Z + acc_array
        
#         if i % 100 == 0:
#             print(f"Processed facility {i}/{len(os.listdir(tt_dir))}, max Z = {np.max(Z)}")
        
#     S = np.zeros_like(mask_arr_cop)
#     U = U/Z

#     for i, file in enumerate(os.listdir(tt_dir)):
#         if not file.endswith('.tif'):
#             continue
#         dist_rast = rasterio.open(tt_dir+file)
#         dist_arr = dist_rast.read()[0]


#         dist_arr = np.nan_to_num(dist_arr, nan = np.inf)

#         # Compute accessibility using vectorized operations
#         p_norm = exp_impedance_fun(dist_arr, beta)
#         p_norm = np.nan_to_num(p_norm, nan = 0)

#         # set places where the friction surface is infinite to nan
#         p_norm[~np.isfinite(mask_arr)] = np.nan
        
#         p_norm = p_norm/Z
#         S -= p_norm*np.nan_to_num(np.log(p_norm), nan = 0)
        
#         if i % 100 == 0:
#             print(f"Processed facility {i}/{len(os.listdir(tt_dir))}")
    
#     S[~np.isfinite(S)] = 0
#     U[~np.isfinite(U)] = 0

#     S[~np.isfinite(mask_arr_cop)] = np.nan
#     U[~np.isfinite(mask_arr_cop)] = np.nan
#     return S, U

def compute_mod_entropy_acc(beta, mask_arr, tt_dir, upper_bound):
    """
    Computes entropy-based accessibility (S) and Average Travel-Time (internal energy) (U) for a given population 
    and travel time rasters. In this version, we introduce an upper bound for the travel time which sets the probability
    of facilities beyond this threshold to zero.

    This function processes travel time raster files in the specified directory, calculates impedance-based 
    accessibility for each facility, and computes two metrics:
    - `S` (Entropy-based accessibility): Measures diversity in accessibility across multiple facilities.
    - `U` (Average Travel Time): Measures the average travel time across all facilities according to travel behaviour.

    Parameters:
    ----------
    beta : float
        The impedance parameter for the exponential decay function, which controls how accessibility 
        decreases with increasing travel time.
    mask_arr : np.ndarray
        Friction surface to use as mask for the country border
    tt_dir : str
        The directory containing `.tif` raster files representing travel time data for different facilities.

    Returns:
    -------
    tuple (np.ndarray, np.ndarray)
        S : np.ndarray
            A 2D numpy array representing entropy-based accessibility for each grid cell.
        U : np.ndarray
            A 2D numpy array representing average travel time for each grid cell.
    """
    
    # Precompute the population mask to avoid redundant operations
    mask_arr_cop = np.where(mask_arr < 0, np.inf, mask_arr.copy())
    
    # Precompute the square of the coordinates to avoid repetitive calculation
    
    Z = np.zeros_like(mask_arr_cop)
    # U = np.zeros_like(mask_arr_cop)

    for i, file in tqdm(enumerate(os.listdir(tt_dir)), desc = "Computing entropy 1/2"):
        if not file.endswith('.tif'):
            continue

        dist_rast = rasterio.open(tt_dir+file)
        dist_arr = dist_rast.read()[0]

        # Compute accessibility using vectorized operations
        acc_array = exp_impedance_fun(dist_arr, beta)
        acc_array = np.nan_to_num(acc_array, nan = 0)
        # set places where the friction surface is infinite to nan
        acc_array[~np.isfinite(mask_arr_cop)] = np.nan
        #dealing with nans
        acc_array[~np.isfinite(acc_array)] = 0
        acc_array[dist_arr > upper_bound] = 0

        dist_arr = np.where(dist_arr == np.inf, 0, dist_arr)

        # U += acc_array*dist_arr
        
        Z = Z + acc_array
        
        # if i % 100 == 0:
        #     print(f"Processed facility {i}/{len(os.listdir(tt_dir))}, max Z = {np.max(Z)}")
        
    S = np.zeros_like(mask_arr_cop)
    # U = U/Z

    for i, file in tqdm(enumerate(os.listdir(tt_dir)), desc = "Computing entropy 2/2"):
        if not file.endswith('.tif'):
            continue
        dist_rast = rasterio.open(tt_dir+file)
        dist_arr = dist_rast.read()[0]


        dist_arr = np.nan_to_num(dist_arr, nan = np.inf)

        # Compute accessibility using vectorized operations
        p_norm = exp_impedance_fun(dist_arr, beta)
        p_norm[dist_arr > upper_bound] = 0

        p_norm = np.nan_to_num(p_norm, nan = 0)

        # set places where the friction surface is infinite to nan
        p_norm[~np.isfinite(mask_arr_cop)] = np.nan

        p_norm = p_norm/Z
        S -= p_norm*np.nan_to_num(np.log(p_norm), nan = 0)
        
        if i % 100 == 0:
            print(f"Processed facility {i}/{len(os.listdir(tt_dir))}")
    
    S[~np.isfinite(S)] = 0
    # U[~np.isfinite(U)] = 0

    S[~np.isfinite(mask_arr_cop)] = np.nan
    # U[~np.isfinite(mask_arr_cop)] = np.nan
    return S
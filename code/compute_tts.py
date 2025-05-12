import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.ops import unary_union
import rasterio
from rasterio.transform import rowcol
from skimage.graph import MCP, MCP_Geometric
import time
from rasterio.transform import from_bounds
from geopy.distance import distance
import os
from shapely.geometry import Polygon, MultiPolygon, Point, mapping
from rasterio.mask import mask
from rasterio.warp import reproject, Resampling, calculate_default_transform


'''
This file contains all the needed functions to compute the 
motorized and walking travel times to the markets using
the motorized and walking Friction Surfaces
'''

# align two rasters with the same crs and resolution
# code from: https://gis.stackexchange.com/questions/404533/aligning-rasters-in-rasterio
def align_extent_raster(raster_reference, tiff_to_be_aligned, output_path):
    """
    Aligns the extent of a raster to match that of a reference raster.

    Parameters:
    ----------
    raster_reference : rasterio.io.DatasetReader
        The reference raster whose extent and profile are used.
    tiff_to_be_aligned : rasterio.io.DatasetReader
        The raster that needs to be aligned to the reference raster.
    output_path : str
        Path where the aligned raster will be saved.

    Returns:
    -------
    None
        The aligned raster is saved to `output_path`.
    """
    profile1 = raster_reference.profile
    bounds1 = raster_reference.bounds
    no_data_1 = raster_reference.nodata
        
    transform_2 = tiff_to_be_aligned.transform

    window = rasterio.windows.from_bounds(*bounds1, transform=transform_2)

    smaller_array = tiff_to_be_aligned.read(window=window, boundless=True) 

    profile1.update(nodata= no_data_1)
    
    # Export it to another image 
    with rasterio.open(output_path, 'w', **profile1) as dst:
        dst.write(smaller_array)

def crop_rast_to_country(rast_path, admin0_gdf, country_num:str, is_pop = False):
    """
    Crops a raster to the boundaries of a country defined by a GeoDataFrame.

    Parameters:
    ----------
    rast_path : str
        Path to the raster file to be cropped.
    admin0_gdf : gpd.GeoDataFrame
        GeoDataFrame containing the country's boundaries.

    Returns:
    -------
    tuple (rasterio.io.DatasetReader, np.ndarray)
        - Cropped raster object.
        - Cropped raster array.
    """
    # Read the Raster
    name, ext = os.path.splitext(rast_path)
    
    cropped_path = name+'_cropped_'+country_num+ext
    
    # Read the raster within the window
    country_rast = rasterio.open(rast_path, 'r')
    crs = country_rast.crs
    '''Cropping the raster to the geometry boundaries'''

    geom = admin0_gdf.to_crs(crs).pop('geometry')
    geom = MultiPolygon(list(geom.apply(lambda x: list(x.geoms) if isinstance(x, MultiPolygon) else x).explode()))

    country_poly = unary_union(geom)
    if isinstance(country_poly, Polygon):
        ext_poly = Polygon(country_poly.exterior)
    else:
        ext_poly = country_poly

    clipped_rast, clipped_trans = mask(country_rast, [mapping(ext_poly)], all_touched = True, crop=True)
    if not is_pop:
        clipped_rast[clipped_rast<0] = np.inf

    cropped = country_rast.meta
    cropped.update({"driver": "GTiff",
                    "height": clipped_rast.shape[1],
                    "width": clipped_rast.shape[2],
                    "transform": clipped_trans})

    with rasterio.open(cropped_path, "w", **cropped) as dest:
        dest.write(clipped_rast)
    country_rast.close() 
    
    country_rast_cropped = rasterio.open(cropped_path, 'r')

    pop_arr = country_rast_cropped.read()[0]
    country_rast_cropped.close()
    return(cropped_path, pop_arr)


'''
Computing Travel-Time rasters

For each market point, we compute one raster where, at each pixel, we have the travel time from it to the market in question.

In order to compute more accurate Travel Times we use the following procedure:

1. Read the FS raster with original geographic CRS
2. Convert the FS to UTM CRS (with constant pixel size)
3. Find average pixel size (avg between vertical, horizontal and diagonal)
4. Adjust the costs from mins/metre to min by multiplying by avg pixel size
5. Compute Shortest path from each market to every pixel
6. Save 1 raster for each market

If we want to compute the distance to the market we'll have to find the length (in steps) of the path and multiply it by the avg pixel size.
'''

def transform_to_utm(input_raster, output_raster, dst_crs):
    """
    Reprojects a raster to the UTM CRS.

    Parameters:
    ----------
    input_raster : str
        Path to the input raster file.
    output_raster : str
        Path to the output raster file in UTM CRS.
    dst_crs : CRS
        Target UTM coordinate reference system.

    Returns:
    -------
    None
        The transformed raster is saved to `output_raster`.
    """
    with rasterio.open(input_raster) as src:
        # Estimate the UTM CRS for the raster
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds
        )
        
        # Update metadata
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        
        # Write the transformed raster
        with rasterio.open(output_raster, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):  # Loop over all raster bands
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest
                )


def calculate_pixel_size_in_meters(transform, crs):
    """
    Calculates the average pixel size in meters based on the raster's affine transformation and CRS.

    Parameters:
    ----------
    transform : Affine
        Affine transformation of the raster.
    crs : CRS
        Coordinate reference system of the raster.

    Returns:
    -------
    float
        Average pixel size in meters.
    """
    # Get the center coordinates of a pixel
    lon1, lat1 = transform * (0.5, 0.5)  # Top-left corner + half pixel shift
    lon2, lat2 = transform * (1.5, 0.5)  # Shift one pixel right
    lon3, lat3 = transform * (0.5, 1.5)  # Shift one pixel down

    # Create GeoDataFrame with the three points
    points_gdf = gpd.GeoDataFrame(
        {'geometry': [Point(lon1, lat1), Point(lon2, lat2), Point(lon3, lat3)]},
        crs=crs
    )
    
    # Estimate UTM CRS and reproject points
    # utm_crs = points_gdf.estimate_utm_crs()
    # points_gdf = points_gdf.to_crs(utm_crs)
    
    # Extract reprojected coordinates
    p1 = points_gdf.geometry.iloc[0]
    p2 = points_gdf.geometry.iloc[1]
    p3 = points_gdf.geometry.iloc[2]
    
    # Calculate Euclidean distances in the projected CRS
    pixel_size_x = p1.distance(p2)  # Horizontal distance in meters
    pixel_size_y = p1.distance(p3)  # Vertical distance in meters

    # Use the average of the two as the pixel size in meters
    return (pixel_size_x + pixel_size_y) / 2
    
# --- Step 3: Generate a raster with travel times to target using MCP.find_costs ---
def generate_travel_time_raster(friction_data, transform, target_coords, pixel_size):
    """
    Generates a raster showing travel times to a target location.

    Parameters:
    ----------
    friction_data : np.ndarray
        Array of friction values (e.g., travel time per unit distance).
    transform : Affine
        Affine transformation for the raster.
    target_coords : tuple
        Geographic coordinates of the target location.
    pixel_size : float
        Average pixel size in meters.

    Returns:
    -------
    np.ndarray
        Raster of travel times to the target location.
    """
    # Get the pixel coordinates for the target destination
    row, col = rowcol(transform, *target_coords)
    
    adjusted_friction_data = friction_data * pixel_size  # Now in minutes per pixel step

    # Create an MCP object, which allows us to find the minimum cost path efficiently
    mcp = MCP_Geometric(adjusted_friction_data)
    target_pixel = [[row, col]]
    # Use the `find_costs` method to calculate the minimum cost (travel time) from the target
    travel_time_raster, traceback = mcp.find_costs(starts = target_pixel)

    # rows, cols = travel_time_raster.shape
    return travel_time_raster
            

# --- Step 4: Visualize the travel time raster ---
def visualize_travel_time_raster(travel_time_raster):
    plt.figure(figsize=(10, 8))
    plt.imshow(travel_time_raster[0], cmap='hot')
    plt.title('Travel Time to Target (Minutes)')
    plt.colorbar(label='Travel Time (Minutes)')
    plt.show()

# reproject from UTM to Geographic
def reproject_to_geographic(src_array, src_transform, src_crs, dst_crs, dst_width, dst_height, dst_transform):
    """
    Reprojects a raster from UTM CRS to geographic CRS.

    Parameters:
    ----------
    src_array : np.ndarray
        Array containing the raster data.
    src_transform : Affine
        Affine transformation of the source raster.
    src_crs : CRS
        Coordinate reference system of the source raster.
    dst_crs : CRS
        Target geographic coordinate reference system.
    dst_width : int
        Width of the destination raster.
    dst_height : int
        Height of the destination raster.
    dst_transform : Affine
        Affine transformation of the destination raster.

    Returns:
    -------
    np.ndarray
        Reprojected raster array.
    """
    # Calculate the bounding box of the input raster in geographic coordinates
    src_bounds = rasterio.transform.array_bounds(
        src_array.shape[0], src_array.shape[1], src_transform
    )
    src_array[np.isinf(src_array)] = np.nan
    # Calculate the destination transform using the target height/width
    # dst_transform = from_bounds(
    #     *src_bounds,  # Use source bounds
    #     width=dst_width,  # Known output width
    #     height=dst_height  # Known output height
    # )
    # Create an empty array for the reprojected data
    dst_array = np.empty((dst_height, dst_width), dtype=src_array.dtype)

    # Reproject the data
    reproject(
        source=src_array,
        destination=dst_array,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        resampling=Resampling.nearest  # Use appropriate resampling
    )
    # print(np.nanmin(dst_array))
    return dst_array

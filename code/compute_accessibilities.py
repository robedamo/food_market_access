import geopandas as gpd
import argparse
import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import rasterio
import scipy
import pickle
from tqdm import tqdm

from shapely.geometry import Polygon, MultiPolygon, Point, mapping

# local files

import centroids 
import compute_tts
import accessibility_metrics
# codes for the countries are taken from openAfrica's dataset (https://open.africa/dataset/africa-shapefiles/resource/dcdadd25-0137-4c93-ae5a-82b39d424d60)

def load_market_data(country_code:str, overlap_rad = 200):
    '''  
    Loads and joins market data from OSM and WFP 

    VARIABLES:
    --------------------------------------------------
    country_code (string): iso3 code of the country to load (lowercase)
    overlap_rad (float): overlap radious in meters

    OUTPUT:

    facility_gdf: Joined GeoDataFrame of OSM and WFP markets
    '''

    # read OSM data
    market_OSM = gpd.read_file('../shared_data/africa_markets/markets/'+country_code+'_markets_shops.geojson')
    print('Number of markets fron OSM:', len(market_OSM))
    crs = market_OSM.crs
    utm_crs = market_OSM.estimate_utm_crs()

    market_OSM = market_OSM.to_crs(utm_crs)
    # get centroids if markets are polygons
    market_OSM['geometry'] = market_OSM['geometry'].centroid
    market_OSM = market_OSM.to_crs(crs)

    # Load WFP data

    # Reading WFP Markets and getting the geometry column out of the .csv file

    market_WFP = pd.read_csv('../shared_data/markets_MFI_africa.csv')
    market_price = pd.read_csv('../shared_data/markets_price_africa.csv')

    market_WFP = pd.concat([market_WFP, market_price]).drop_duplicates(subset = 'MarketId', keep = 'first')

    market_WFP['geometry'] = [Point(xy) for xy in zip(market_WFP.Longitude, market_WFP.Latitude)]
    
    market_WFP = gpd.GeoDataFrame(market_WFP, geometry = 'geometry')
    
    market_WFP.drop(['Latitude', 'Longitude'], inplace = True, axis = 1)
    market_WFP = market_WFP.set_crs(crs)
    market_WFP = market_WFP.to_crs(utm_crs)

    print(len(market_WFP))

    # load the country borders
    country_borders = gpd.read_file('../shared_data/africa_markets/borders/'+country_code+'.geojson')
    country_borders = country_borders.to_crs(utm_crs)

    # Get the WFP markets from the country considered

    market_WFP_country= market_WFP.sjoin(country_borders)
    market_WFP_country.reset_index(inplace = True)
    
    print(f'Number of markets fron WFP in {country_borders.iloc[0].ADM0_NAME}:', len(market_WFP_country))

    # Now we merge the WFP and OSM datasets and set an overlap radious of 200m
    # go to UTM crs in both gdfs

    market_OSM = market_OSM.to_crs(utm_crs)
    OSM_points_x = np.array([point.x for point in market_OSM['geometry']])
    OSM_points_y = np.array([point.y for point in market_OSM['geometry']])
    same_points = []
    
    for i in range(len(market_WFP_country)):
        new_market = market_WFP_country['geometry'].iloc[i]
        dists_x = OSM_points_x - new_market.x
        dists_y = OSM_points_y - new_market.y

        dists = np.sqrt(dists_x**2 + dists_y**2)

        # if the dist. between markets is less than ri they are the same
        if np.min(dists) < overlap_rad:
            same_points.append(i)
            
    # drop the repeated points
    market_WFP_dropped = market_WFP_country.drop(index = same_points, axis = 0)
    
    # set a market ID
    # MFI_WFP_dropped['ID'] = list(range(max_ID+1, max_ID+len(MFI_WFP_dropped)+1))
    
    # concatenate the WFP and OSM dataframes
    market_gdf =  pd.concat([market_OSM,market_WFP_dropped])
    market_gdf = market_gdf.to_crs(crs)
    market_gdf['ID'] = np.arange(0, len(market_gdf))
    market_gdf = market_gdf[['ID', 'geometry']]
    
    facility_gdf = market_gdf.to_crs(utm_crs).sjoin(country_borders.to_crs(utm_crs), how = 'left')

    facility_gdf = facility_gdf.to_crs(crs)
    facility_gdf = facility_gdf.reset_index()
    facility_gdf = facility_gdf.drop_duplicates(subset = 'ID')
   
    facility_gdf = facility_gdf[['ID', 'geometry']]

    
    print('Total number of points', len(facility_gdf))
    print('Number of overlapping points at less than {0}m:'.format(overlap_rad), len(same_points))

    return(facility_gdf, country_borders.to_crs(crs))


# Adding colorbars
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if x!=0:
        return r'${} \times 10^{{{}}}$'.format(a, b)
    else:
        return '0'
    
### MAIN PROGRAM ###

def main(country_code, threshold, tmode):
    if tmode != 'moto' and tmode != "walk":
        msg = f"Chose a valid transport mode. Options walk for walking and moto for motorized. {tmode} is not valid"
        raise Exception(msg)
    try:
        threshold = float(threshold)
    except:
        print("Provide the threshold as a number in string format (e.g. '60')")
        return

    if len(country_code) != 3:
        msg = f"Input the country code in a valid iso3 format (e.g. 'eth' for Ethiopia)"
        raise Exception(msg)
    
    # change to lowercase format to prevent errors

    country_code = country_code.lower()
    if len(country_code) != 3:
        msg = f"Input the country code in a valid iso3 format (e.g. 'eth' for Ethiopia)"
        raise Exception(msg)
    

    facility_gdf, country_gdf = load_market_data(country_code)

    utm_crs = facility_gdf.estimate_utm_crs()
    if not os.path.exists('./computed_centroids/'+country_code+'_centroids.geojson'):    

        # Spatial join: Finds points within the polygon
        joined_gdf = gpd.sjoin(facility_gdf.to_crs(utm_crs), country_gdf.to_crs(utm_crs), predicate="within")
        points = np.array([[point.x,point.y]  for point in joined_gdf.to_crs(utm_crs).geometry])
            # if len(points) == 0:
            #     empty_polys+=1
            #     print('Empty polygon', empty_polys)
            #     continue
        # overlap radious (15 min at walking pace of 1.4 m/s)
        ri = 1.4*15*60 
        if len(points) == 1:
            centers = points
            clusters = np.array([[points[0][0],points[0][1],0]])
        for i in range(2,len(points)):
            if not i%10:
                print(i)
            valid, centers, clusters = centroids.validate_solution(ri, *centroids.create_clusters(i,points))
            if valid:
                print(i)
                break
        
        centroids_gdf, facility_gdf = centroids.create_centroids_gdf(facility_gdf, centers, clusters)

        centroids_gdf.to_file('./computed_centroids/'+country_code+'_centroids.geojson')
        print('Centroids computed:', str(len(centroids_gdf))+' centroids')
    else:
        centroids_gdf = gpd.read_file('./computed_centroids/'+country_code+'_centroids.geojson')

# If needed you can use the worldpop data
    # wpp_path = '../shared_data/worldpop_rasters/'+country_code

    # for file in os.listdir(wpp_path):
    #     if file.endswith("UNadj.tif"):
    #         wpp_path+='/'+file

    # cropped_path_wpp, pop_arr = compute_tts.crop_rast_to_country(wpp_path, country_gdf, country_num=country_code, is_pop = True)

    # Checking that the raster and the points match...

    # Loading Friction Surface
    print('preparing friction surface to compute travel times ...')
    if tmode == 'moto':
        fs_path = '../shared_data/friction_surfaces/2020_motorized_friction_surface.geotiff'
        suffix = '_moto'
    if tmode == 'walk':
        fs_path = '../shared_data/friction_surfaces/2020_walking_only_friction_surface.geotiff'
        suffix = '_walk'
    try:
        cropped_path, fs_arr = compute_tts.crop_rast_to_country(fs_path, country_gdf, country_num=country_code, is_pop = False)
    except Exception as e:
        raise e
    # transform FS to UTM to compute TT rasters

    name, ext = os.path.splitext(fs_path)
    out_path_utm = name+'_'+country_code+'_utm'+ext
    
    compute_tts.transform_to_utm(cropped_path, out_path_utm, utm_crs)

    # File path to the friction raster
    root_dest = './computed_tts/'
    dest_path = root_dest+country_code
        # Compute travel times if they have not been computed for this country
    if not os.path.exists(root_dest): 
        os.mkdir(root_dest)

    if not os.path.exists(dest_path): 
        os.mkdir(dest_path)
        # Load the raster
        with rasterio.open(out_path_utm) as src:
            friction_data = src.read(1)  # Read the first band (assuming friction values)
            friction_data[friction_data<0] = np.inf
            src_transform = src.transform  # Affine transform for the raster
            # width, height = src.width, src.height  # Dimensions of the raster
            src_crs = src.crs
            # metadata = src.meta

        # Calculate pixel size in meters for geographic CRS
        pixel_size = compute_tts.calculate_pixel_size_in_meters(src_transform, src_crs)


        # Target dimensions (known height and width)
        target_width = fs_arr.shape[1]  
        target_height = fs_arr.shape[0] 
        fs_src = rasterio.open(cropped_path)
        dst_transform = fs_src.transform
        dst_crs = fs_src.crs
        fs_src.close()
        for i in tqdm(range(len(centroids_gdf)), desc = "Computing Travel Time Rasters"):
            
            # if i%100==0:
            #     print(i)
            target = centroids_gdf['geometry'].to_crs(src_crs).iloc[i]
            target_coords = [target.x, target.y]
            ID = centroids_gdf.ID.iloc[i]
            # Generate the travel time raster using MCP.find_costs
            travel_time_raster = compute_tts.generate_travel_time_raster(friction_data, src_transform, target_coords, pixel_size)

            # Transform, and CRS in UTM
            reprojected_array = compute_tts.reproject_to_geographic(travel_time_raster, src_transform, 
                                                        src_crs, dst_crs, target_width, 
                                                        target_height, dst_transform)
            

            # Save the result to a GeoTIFF
            with rasterio.open(dest_path+'/'+str(ID)+'.tif',
                "w",
                driver="GTiff",
                height=target_height,
                width=target_width,
                count=1,
                dtype=reprojected_array.dtype,
                crs=dst_crs,
                transform=dst_transform) as dst:
            
                dst.write(reprojected_array, 1)

    else:
        print(f'Travel time rasters already computed for country {country_code}')

    '''COMPUTING ACCESSIBIITIES'''

    R = 60 # min
    beta = -np.log(0.01)/R # accessibility = 0.01 when r = 60 min
    tt_dir = dest_path+'/'

    # print('Computing Gravity Acc.')
    # access_mat = accessibility_metrics.compute_grav_acc(beta, pop_arr, tt_dir)
    
    # print('Computing Cumulative Acc.')
    # access_mat_cumul =  accessibility_metrics.compute_grav_acc_cumul(threshold, fs_arr, tt_dir)

    print('Computing Pure Cumulative Acc.')
    access_mat_pure_cumul =  accessibility_metrics.compute_cumulative(threshold, fs_arr, tt_dir)

    # print('Computing Entropy')
    # entropy, mean_dist = accessibility_metrics.compute_entropy_acc(beta, fs_arr, tt_dir)

    print('Computing Mod Entropy with 4h upper bound')
    entropy_mod = accessibility_metrics.compute_mod_entropy_acc(beta, fs_arr, tt_dir, 4*60)

    print('Computing Shortest TT.')
    shortest_tt_rast = accessibility_metrics.shortest_travel_time(tt_dir)

    # print('max tt', np.nanmax(shortest_tt_rast))
    '''SAVING ACCESSIBIITIES'''

    if not os.path.exists('./computed_acc/'+country_code): 
        os.makedirs('./computed_acc/'+country_code)


    mod_en_path = './computed_acc/'+country_code+f'/entropy_mod_4h{suffix}.tif'


    shortest_tt_path = './computed_acc/'+country_code+f'/shortest{suffix}_tt.tif'
    acc_path_pure_cumul = './computed_acc/'+country_code+f'/acc_pure_cumul{suffix}.tif'
    
    fs_src = rasterio.open(cropped_path)
    dst_transform = fs_src.transform
    dst_crs = fs_src.crs
    fs_src.close()
    # Write the entropy 
    
    with rasterio.open(mod_en_path,"w",driver="GTiff", height=entropy_mod.shape[0],
                       width=entropy_mod.shape[1], count=1, dtype=entropy_mod.dtype,
                       crs=dst_crs, transform=dst_transform) as dst:
        
        dst.write_band(1, entropy_mod) 

    # Write the travel time 
    
    with rasterio.open(shortest_tt_path,"w",driver="GTiff", height=shortest_tt_rast.shape[0],
                       width=shortest_tt_rast.shape[1], count=1, dtype=shortest_tt_rast.dtype,
                       crs=dst_crs, transform=dst_transform) as dst:
        
        dst.write_band(1, shortest_tt_rast) 

    # Write the cumulative 

    with rasterio.open(acc_path_pure_cumul,"w",driver="GTiff", height=access_mat_pure_cumul.shape[0],
                       width=access_mat_pure_cumul.shape[1], count=1, dtype=access_mat_pure_cumul.dtype,
                       crs=dst_crs, transform=dst_transform) as dst:

        dst.write_band(1, access_mat_pure_cumul) 
        
    # finally delete the travel time rasters due to memory concerns
    #remove the friction raster utm file to avoid trash
    os.remove(out_path_utm)
    os.remove(cropped_path)
    shutil.rmtree(tt_dir)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process country data")
    parser.add_argument("country_code", type=str, help="Number of the country to process")
    parser.add_argument("threshold", type=str, help="threshold travel time for pure cumulative accessibility")
    parser.add_argument("tmode", type=str, help="Transport mode ('moto' or 'walk')")

    args = parser.parse_args()

    main(args.country_code, args.threshold, args.tmode)

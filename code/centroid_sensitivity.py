import geopandas as gpd
import argparse
import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import rasterio
import scipy
from tqdm import tqdm
import csv
from shapely.geometry import Polygon, MultiPolygon, Point, mapping

# local files

import centroids 

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


def main(country_code:str):

    radius_ls = [5, 10, 15, 20, 25, 30, 35, 40]
    lens_clusters = []
    facility_gdf, country_gdf = load_market_data(country_code)
    utm_crs = facility_gdf.estimate_utm_crs()

    # Spatial join: Finds points within the polygon
    joined_gdf = gpd.sjoin(facility_gdf.to_crs(utm_crs), country_gdf.to_crs(utm_crs), predicate="within")
    points = np.array([[point.x,point.y]  for point in joined_gdf.to_crs(utm_crs).geometry])
    lens_clusters.append(len(points))
    for radius in tqdm(radius_ls, desc='Radius', unit='m'):
        print(f'Radius: {radius}min')
        # overlap radius (15 min at walking pace of 1.4 m/s)
        ri = 1.4*radius*60 
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
        
        centroids_gdf, _ = centroids.create_centroids_gdf(facility_gdf.copy(), centers, clusters)
        lens_clusters.append(len(centroids_gdf))

    with open(f'./cluster_sensitivity/{country_code}_clusters.csv', 'w') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow([0, lens_clusters[0]])
        for i in range(len(radius_ls)):
            csv_writer.writerow([radius_ls[i], lens_clusters[i]])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process country data")
    parser.add_argument("country_code", type=str, help="Number of the country to process")
    print('in main')
    args = parser.parse_args()

    main(args.country_code)
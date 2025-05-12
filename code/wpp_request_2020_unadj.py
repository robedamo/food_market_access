#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 11:17:59 2024

@author: benassi
"""

import geopandas as gpd
import os
import requests
import rasterio
from rasterio.merge import merge, copy_sum
import numpy as np
from bs4 import BeautifulSoup
from urllib.parse import urljoin

gdf = gpd.read_file('/Users/benassi/Documents/ISI/wfp-micronutrient-2024/shared_data/africa_borders/afr_g2014_2013_0.shp')



def download_worldpop_data(country_code, year, output_dir="worldpop_data"):
    """
    Download WorldPop UN-adjusted population data for a specific country and year.

    Args:
        country_code (str): ISO3 country code (e.g., "KEN" for Kenya).
        year (int): Year of population data (e.g., 2020).
        output_dir (str): Directory to save the downloaded file.

    Returns:
        str: Path to the downloaded file.
    """
    lower_countr_code = country_code.lower()

    base_url = "https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj"
    file_name = f"{lower_countr_code}_ppp_{year}_1km_Aggregated_UNadj.tif"
    download_url = f"{base_url}/{year}/{country_code}/{file_name}"
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, file_name)

    try:
        print(f"Downloading data from {download_url}...")
        response = requests.get(download_url, stream=True)
        response.raise_for_status()  # Raise HTTPError for bad responses

        # Save the file
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
        print(f"Data downloaded successfully: {output_path}")
        return output_path

    except requests.exceptions.RequestException as e:
        print(f"Error downloading data: {e}")
        return None

def download_worldpop_demographics_data(country_code, year, output_dir="worldpop_data"):
    """
    Download WorldPop UN-adjusted population data for a specific country and year.

    Args:
        country_code (str): ISO3 country code (e.g., "KEN" for Kenya).
        year (int): Year of population data (e.g., 2020).
        output_dir (str): Directory to save the downloaded file.

    Returns:
        str: Path to the downloaded file.
    """

    base_url = f"https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_1km/unconstrained/{year}/{country_code}/"

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    # dir_ls = os.listdir(output_dir)
    # if np.any([file.startswith(country_code+"_merged_raster_male_") for file in dir_ls]):
    #     print(f'{country_code} already computed.')
    #     return
    # Get the HTML content of the page
    response = requests.get(base_url)
    if response.status_code != 200:
        print("Failed to retrieve the webpage")
        exit()

    # Parse the page
    soup = BeautifulSoup(response.text, "html.parser")

    # Find all links in the page
    links = soup.find_all("a")

    # Download files
    for link in links:
        file_url = urljoin(base_url, link.get("href"))  # Get absolute URL

        if file_url.endswith((".tif")):  # Adjust file extensions if needed
            file_name = os.path.join(output_dir, os.path.basename(file_url))
            
            print(f"Downloading {file_url}...")
            # if os.path.exists(file_name):
            #     continue
            # Download the file
            with requests.get(file_url, stream=True) as file_response:
                if file_response.status_code == 200:
                    with open(file_name, "wb") as f:
                        for chunk in file_response.iter_content(chunk_size=1024):
                            f.write(chunk)
                    print(f"Saved: {file_name}")

                else:
                    print(f"Failed to download {file_url}")

    sum_rasters(output_dir, country_code)
    print(f'Completed {country_code}')
    
def sum_rasters(directory:str, country_code:str):

    # find all male files
    male_files = [f for f in os.listdir(directory) if f[4] =='m']
    female_files = [f for f in os.listdir(directory) if f[4] =='f']

    # group into age categories
    age_ls =[10, 20, 30, 40, 50, 60, 70]

    for i, age in enumerate(age_ls):
        if i == 0:
            age_l = 0
        else:
            age_l = age_ls[i-1]+5

        male_age_gap = []

        for file in male_files:
            try:
                if file[6:8].isnumeric() and int(file[6:8]) > age_l and int(file[6:8]) <= age:
                    male_age_gap.append(file)
                else:
                    continue
            except (ValueError, IndexError):
                print(file[6])
                if int(file[6]) > age_l and int(file[6]) <= age:
                    male_age_gap.append(file)
                else:
                    continue

        female_age_gap = []
        
        for file in female_files:
            try:
                if file[6:8].isnumeric() and int(file[6:8]) > age_l and int(file[6:8]) <= age:
                    female_age_gap.append(file)
                else:
                    continue
            except (ValueError, IndexError):
                if int(file[6]) > age_l and int(file[6]) <= age:
                    female_age_gap.append(file)
                else:
                    continue
    
        # Open all rasters

        src_male_files = [rasterio.open(os.path.join(directory,raster)) for raster in male_age_gap]
        merged_raster_male, _ = merge(src_male_files, method = 'sum')

        src_female_files = [rasterio.open(os.path.join(directory,raster)) for raster in female_age_gap]
        merged_raster_female, _ = merge(src_female_files, method = 'sum')

        # Get metadata from one of the source files
        src_meta = src_male_files[0].meta.copy()
        # Update metadata for the merged raster
        src_meta.update({
            "height": merged_raster_male.shape[1],
            "width": merged_raster_male.shape[2],
            "count": merged_raster_male.shape[0]
        })
        
        # Save the merged raster
        output_path = os.path.join(directory, country_code+f"_merged_raster_male_{age_l}_{age+5}.tif")
        # if not os.path.exists(output_path):
        with rasterio.open(output_path, "w", **src_meta) as dst:
            dst.write(merged_raster_male)
        
        # Save the merged raster
        output_path = os.path.join(directory, country_code+f"_merged_raster_female_{age_l}_{age+5}.tif")
        # if not os.path.exists(output_path):
        with rasterio.open(output_path, "w", **src_meta) as dst:
            dst.write(merged_raster_female)
        
    # delete original files
    for file in male_files:
        os.remove(os.path.join(directory,file))

    for file in female_files:
        os.remove(os.path.join(directory,file))

# no_iso3 = []
# for row in range(len(gdf)):
#     iso3 = gdf.iloc[row]['ISO3']
#     num_code = gdf.iloc[row]['ADM0_CODE']
#     if iso3 == None:
#         country_nm = gdf.iloc[row]['ADM0_NAME']
#         no_iso3.append(country_nm)
#         continue
#     out_dir = f'/Users/benassi/Documents/ISI/wfp-micronutrient-2024/shared_data/worldpop_rasters/demographics/{num_code}/'
#     download_worldpop_demographics_data(iso3, 2020, output_dir=out_dir)

out_dir = f'/Users/benassi/Documents/ISI/wfp-micronutrient-2024/shared_data/worldpop_rasters/demographics/{268}/'
print(gdf[gdf['ADM0_CODE'] == 268]['ADM0_NAME'])
download_worldpop_demographics_data('ESH', 2020, output_dir=out_dir)
        
# print('countries without iso3:', no_iso3)

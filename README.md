# A large-scale analysis of food market spatial accessibility in Africa

## Table of Contents
1. [Citing](#citing)
2. [Abstract](#abstract)
3. [Description](#description)
4. [Usage](#usage)
5. [Authors and acknowledgment](#authors-and-acknowledgment)
6. [Requirements](#requirements)
7. [License](#license)

## 1. Citing

## 2. Abstract

Access to food is a fundamental human right and is comprised in the United Nations Sustainable Development Goals. However, market access remains limited, particularly in low- and middle-income countries, where research on market accessibility is scarce and often relies on small-scale surveys with limited quantitative scope. In this study, we present a large-scale analysis of spatial market accessibility across Africa, leveraging open data from OpenStreetMap and the United Nations World Food Programme. We evaluate accessibility through three key dimensions: travel time (walking or motorized), market availability, and spatial distribution. Our findings reveal consistent patterns across all three metrics, highlighting stark disparities between urban and rural areas, with the latter, often less wealthy regions, experiencing severe access limitations. Furthermore, we examine the link between market access and food security, identifying a significant yet moderate correlation. These results underscore the critical role of spatial market access in food security assessments. Our scalable approach provides a valuable framework for identifying regions with the greatest accessibility challenges, serving as a foundation for future targeted research and policy interventions.

## 3. Description

This repository contains the code needed to compute all the spatial accessibility metrics for any African country as well as the script that produces the figures present in the published manuscript. We use three spatial accessibility definitions:

- **Cumulative accessibility**: Number of markets within $T$ minutes of travel time.
- **Travel time to closest market**: Travel time computed from every pixel to the closest market both in motorized and walking transport modes.
- **Entropy**: Entropy of a poximity dependent travel time distribution (see the manuscript for more details).

The main script computing all the accessibility metrics is `compute_accessibilities.py` inside the `code` folder. This file calls three libraries:

1. `compute_tts.py`: Computes the travel times from each pixel to every market point using a friction surface (motorized or walking). The main funciton produces one raster for each market, with the time from each pixel to that point. After the main script finishes, it deletes the rasters to avoid accumulating heavy files.
2. `centroids.py`: Takes a set of locations (food shops and markets) and groups them them in clusters of a maximum of 15 min. at a walking velocity.
3. `accessibility_metrics.py`: Contains all the accessibility metrics used (travel time, cumulative and entropy) as well as some extras (pure entropy, shen and hansen type accessibilities.).

The main function from `compute_accessibilities.py` takes in 3 variables. First the country code (numeric), secondly the cumulative accessibility travel time threshold (in minutes) and thirdly the travel mode (motorized or walking). The code saves the accessibility matrices in `.pickle` format in `.\code\computed_acc\{country number}` for the specified country.

This project uses data from [OpenStreetMap (OSM)](https://www.geofabrik.de/geofabrik/openstreetmap.html) and the [World food Programme (WFP)](https://www.wfp.org/). We also make use of the friction surfaces from the [Malaria Atlas Project](https://malariaatlas.org/project-resources/accessibility-to-healthcare/). Check the [Data and Extraction](#data-and-extraction)senction for more information on the data. You can find all the needed data to compute the plots of the paper in the `shared_data` folder.

## 4. Data and Extraction

### 4.1 OSM

To extract the OSM data, we used [Geofabrik](https://download.geofabrik.de/)'s downoad server. In our case, we downoaded the african dataset in the `osm.pbf` format. You can find the bash extraction script in `shared_data/extract_countries.sh`. You will need the [`osmium`](https://osmcode.org/osmium-tool/) tool to extract the markets from the geofabrik file using the provided script. You will also need [`ogr2ogr`](https://gdal.org/en/stable/programs/ogr2ogr.html) from [gdal](https://gdal.org/en/stable/).

We extract the following amenities:

- [`marketplace`](https://wiki.openstreetmap.org/wiki/Tag:amenity%3Dmarketplace)
- [`supermarket`](https://wiki.openstreetmap.org/wiki/Tag:shop%3Dsupermarket)
- [`bakery`](https://wiki.openstreetmap.org/wiki/Tag:shop%3Dbakery)
- [`butcher`](https://wiki.openstreetmap.org/wiki/Tag:shop%3Dbutcher)
- [`convenience`](https://wiki.openstreetmap.org/wiki/Tag:shop%3Dconvenience)
- [`dairy`](https://wiki.openstreetmap.org/wiki/Item:Q6186)
- [`food`](https://wiki.openstreetmap.org/wiki/Tag:shop%3Dfood)
- [`farm`](https://wiki.openstreetmap.org/wiki/Tag:shop%3Dfarm)

### 4.2 WFP

We use the WFP's [Prices](https://dataviz.vam.wfp.org/economic/prices?current_page=1&theme=10) and [Market assessment](https://dataviz.vam.wfp.org/economic/market-assessment?current_page=1&theme=30) data, downloaded from their sites for all african countries. You can find the data in `shared_data\markets_MFI_africa.csv` and `shared_data\markets_price_africa.csv`.

### 4.3 Worldpop

For analysis purposes we use the [WorldPop UN-adjusted population counts at 1km resolution](https://hub.worldpop.org/geodata/listing?id=75) (30-arc seconds).

### 4.4 Friction Surface
We use the friction surfaces from the [Malaria Atlas Project](https://malariaatlas.org/project-resources/accessibility-to-healthcare/) to compute the travel times to markets. These surfaces are based on topographic data as well as transportation networks. The units are seconds to travel 1 meter.

### 4.5 Rural-urban classification

For the analysis on the relation between rurality and market access we use a 30-arcsec raster from the [Urban-Rural Catchment Areas (URCA)](https://figshare.com/articles/dataset/Urban-rural_continuum/12579572) Project, with a 30-level classificaiton.

## 5. Usage

To call `compute_accessibilities.py`, you need to input three variables:
- **country_code**: A string with the iso3 code of the country (e.g. `'eth'` for Ethiopia)
- **threshold**: The travel time threshold (in minutes) for the cumulative accessibility in `str` format (e.g. `'30'`for 30 minutes)
- **tmode**: transportation mode. Use `'moto'`for motorized or `'walk'`for walking-only.

e.g. `python compute_accessibilities.py 'eth' '30' 'moto'`

Note that to run the code, four data files are needed:
1. The market locations in `shared_data/africa_markets/markets/{iso3}_markets_shops.geojson` from OSM.
2. The Market locations from WFP. For Africa `shared_data/markets_MFI_africa.csv` and `/shared_data/markets_price_africa.csv`.
3. The borders of the country in `shared_data/africa_markets/borders/{iso3}.geojson`
4. The Friction Surfaces (FS). In this case the global FS are in `shared_data/friction_surfaces`

The data in this repository allows computing the metrics for all African countries. If you need to use the code for non-African countries download the required data (market locations from OSM and WFP and country borders). Make sure to put the new data in the right directories.

The file in `code/manuscript_plots.ipynb` reproduces all plots present in the paper.

## 6. Requirements
```
beautifulsoup4==4.13.4
geopandas==0.14.0
geopy==2.4.0
matplotlib==3.7.1
numpy==1.24.3
pandas==1.5.3
rasterio==1.3.9
Requests==2.32.3
scikit_learn==1.3.0
scipy==1.10.1
Shapely==2.1.0
skimage==0.0
tqdm==4.65.0
```

## 7. License
This project is licensed under the [MIT License](https://mit-license.org/).

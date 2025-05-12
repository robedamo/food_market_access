# TO EXTRACT ONLY ONE COUNTRY

Define paths
# Put your pbf file here
AFRICA_PBF="../data/africa-latest.osm.pbf"
# Put your geodataframe with the country borders here
SHAPEFILE="./africa_borders/afr_g2014_2013_0.shp"

# Define the country code (in this case numeric coode for South Africa)
COUNTRY="160"

# Loop through each country and process
echo "Processing ${COUNTRY}"
# Export the country's boundary as GeoJSON
ogr2ogr -f "GeoJSON" "./africa_markets/borders/${COUNTRY}.geojson" $SHAPEFILE -where "ADM0_CODE = '${COUNTRY}'"

# Extract the OSM data for the country
osmium extract --polygon="./africa_markets/borders/${COUNTRY}.geojson" $AFRICA_PBF -o "${COUNTRY}.osm.pbf" --overwrite

# Apply tag filters and directly export to GeoJSON
osmium tags-filter -o "${COUNTRY}_markets_shops.osm.pbf" "${COUNTRY}.osm.pbf" \
    nwr/amenity=marketplace \
    nwr/shop=supermarket \
    nwr/shop=bakery \
    nwr/shop=butcher \
    nwr/shop=convenience \
    nwr/shop=dairy \
    nwr/shop=food \
    nwr/shop=farm 

osmium export "${COUNTRY}_markets_shops.osm.pbf" -o "./africa_markets/markets/${COUNTRY}_markets_shops.geojson" --overwrite

rm "${COUNTRY}.osm.pbf"

# Check if output file was created successfully
if [ -s "./africa_markets/markets/${COUNTRY}_markets_shops.geojson" ]; then
    rm "${COUNTRY}_markets_shops.osm.pbf"
    echo "Finished processing ${COUNTRY}"
else
    echo "No data for ${COUNTRY}. Skipping."
    continue
fi

# TO EXTRACT ALL COUNTRIES

# #!/bin/bash

# # Define paths
# AFRICA_PBF="africa-latest.osm.pbf"
# SHAPEFILE="./africa_borders/afr_g2014_2013_0.shp"

# # Get a list of all country names from the shapefile
# COUNTRY_CODES=$(ogrinfo -q -sql "SELECT ADM0_CODE FROM afr_g2014_2013_0" $SHAPEFILE | grep ADM0_CODE | awk '{print $4}')

# mkdir -p "africa_markets/borders"
# mkdir -p "africa_markets/markets"

# # Loop through each country and process
# for COUNTRY in $COUNTRY_CODES; do
#     echo "Processing ${COUNTRY}"
#     # Export the country's boundary as GeoJSON
#     ogr2ogr -f "GeoJSON" "africa_markets/borders/${COUNTRY}.geojson" $SHAPEFILE -where "ADM0_CODE = '${COUNTRY}'"
    
#     # Extract the OSM data for the country
#     osmium extract --polygon="africa_markets/borders/${COUNTRY}.geojson" $AFRICA_PBF -o "${COUNTRY}.osm.pbf"
    
#     # Apply tag filters and directly export to GeoJSON
#     osmium tags-filter -o "${COUNTRY}_markets_shops.osm.pbf" "${COUNTRY}.osm.pbf" \
#     nwr/amenity=marketplace \
#     nwr/shop=supermarket \
#     nwr/shop=bakery \
#     nwr/shop=butcher \
#     nwr/shop=convenience \
#     nwr/shop=dairy \
#     nwr/shop=food \
#     nwr/shop=farm 
   
#     osmium export "${COUNTRY}_markets_shops.osm.pbf" -o "africa_markets/markets/${COUNTRY}_markets_shops.geojson"

#     rm "${COUNTRY}.osm.pbf"

#     # Check if output file was created successfully
#     if [ -s "africa_markets/markets/${COUNTRY}_markets_shops.geojson" ]; then
#         rm "${COUNTRY}_markets_shops.osm.pbf"
#         echo "Finished processing ${COUNTRY}"
#     else
#         echo "No data for ${COUNTRY}. Skipping."
#         continue
#     fi
# done

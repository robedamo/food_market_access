# TO EXTRACT ONLY ONE COUNTRY
COUNTRY='lka'

# Extract the country borders from OSM. Admin_level 2 is country level for OSM.
#osmium tags-filter "${COUNTRY}.osm.pbf" r/boundary=administrative r/admin_level=2 -o "${COUNTRY}_borders.osm.pbf"

#osmium export "${COUNTRY}_borders.osm.pbf" -o "./africa_markets/borders/${COUNTRY}.geojson"

# Check if output file was created successfully
#if [ -s "./africa_markets/borders/${COUNTRY}.geojson" ]; then
#    rm "${COUNTRY}_borders.osm.pbf"
#    echo "Borders extracted for ${COUNTRY}"
#else
#    echo "No borders for ${COUNTRY}. Skipping."
#    continue
#fi

# Extract markets and food shops
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

#rm "${COUNTRY}.osm.pbf"

# Check if output file was created successfully
if [ -s "./africa_markets/markets/${COUNTRY}_markets_shops.geojson" ]; then
    rm "${COUNTRY}_markets_shops.osm.pbf"
    echo "Finished processing ${COUNTRY}"
else
    echo "No data for ${COUNTRY}. Skipping."
    continue
fi


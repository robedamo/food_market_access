#!/bin/bash

# Set the path to the parent directory containing the numbered folders
PARENT_DIR="../shared_data/worldpop_rasters"
CHECK_DIR="./computed_acc"

# Loop through each folder in the parent directory
for folder in "$PARENT_DIR"/*; do
    if [ -d "$folder" ]; then # && [ ! -e "$CHECK_DIR/$(basename "$folder")" ]; then

    	folder_name=$(basename "$folder")
        if [ ${#folder_name} -gt 3 ]; then
            continue
        fi
    # Launch a screen session in detached mode
    	screen_name="screen_$folder_name"
    	echo "Launching screen: $screen_name for folder $folder_name"
    
        # Run the Python script inside the screen
        screen -dmS "$screen_name" bash -c "
                python centroid_sensitivity.py \"$folder_name\""

            #python compute_accessibilities.py \"$folder_name\" "60" "walk""
    fi
done

echo "All screens launched. Use 'screen -ls' to see active sessions."


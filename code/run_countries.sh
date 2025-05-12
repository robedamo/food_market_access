#!/bin/bash

# List of country names

#{6: 'Sudan', 8: 'Angola', 29: 'Benin', 35: 'Botswana', 42: 'Burkina Faso', 45: 'Cameroon', 47: 'Cape Verde', 49: 'Central African Republic', 
#50: 'Chad', 58: 'Comoros', 59: 'Congo', 66: "Côte d'Ivoire", 68: 'Democratic Republic of the Congo', 70: 'Djibouti', 74: 'South Sudan', 
#76: 'Equatorial Guinea', 77: 'Eritrea', 79: 'Ethiopia', 89: 'Gabon', 90: 'Gambia', 94: 'Ghana', 105: 'Guinea-Bissau', 106: 'Guinea', 133: 'Kenya', 
#142: 'Lesotho', 144: 'Liberia', 150: 'Madagascar', 152: 'Malawi', 155: 'Mali', 159: 'Mauritania', 160: 'Mauritius', 170: 'Mozambique', 172: 'Namibia',
#181: 'Niger', 182: 'Nigeria', 217: 'Senegal', 220: 'Seychelles', 221: 'Sierra Leone', 226: 'Somalia', 227: 'South Africa', 235: 'Swaziland', 243: 'Togo', 
#253: 'Uganda', 257: 'United Republic of Tanzania', 270: 'Zambia', 271: 'Zimbabwe', 4: 'Algeria', 43: 'Burundi', 102: 'Abyei', 145: 'Libya', 169: 'Morocco', 
#205: 'Rwanda', 214: 'Sao Tome and Principe', 248: 'Tunisia', 268: 'Western Sahara', 40760: "Hala'ib triangle", 40762: "Ma'tan al-Sarra", 40765: 'Egypt', 
#61013: 'Ilemi triangle'}

#countries=("150" "172" "8")
#countries=("226" "94" "105")
#countries=("160" "76" "170")
# countries=("106" "257" "227")


# Loop through each country and run the Python script
# for country in "${countries[@]}"; do
#     echo "Processing country: $country"
#     python compute_accessibilities.py "$country"
# done

# echo "Finished processing countries: ${countries[*]}"

# Generate a sequence of floats using seq and awk
LC_NUMERIC=C
countries=("235" "74" "77" "205" "270" "90" "42" "70" "214" "271" "35" "50" "29" "47" "142" "58" "181" "159" "268" "145" "49")

# Loop through each travel time in the array
for tt in "${countries[@]}"; do
    echo "Processing travel time: $tt"
    python compute_accessibilities.py "$tt" "30"
done

# Final message
echo "Finished processing travel times: ${arr[*]}"

# Omicron countries split datasets

This script pulls the history of case and sequence counts of datasets that includes the following variant categories: other, Delta, Omicron 21K and Omicron 21L and only includes countries with >100 Omicron 21L sequences that are available at the blab/rt-from-frequency-dynamics repository
 
## Downloading datasets requires cloning rt-from-frequency-dynamics repository

This script outputs case and sequences counts organized by date-stamped folders and can be run using the following script using inputs of location of repository and location of folder where files can be saved

 ```python3 generatre_urls_data.py --location "LOCATION OF REPO" --output "OUTPUT LOCATION"```

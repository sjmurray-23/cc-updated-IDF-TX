"""
Date: 12/4/2025
Time: 11:28 AM
AIM: Calculate the Coefficient of Variation of the Change Factor estimates prior to the interpolation of the results.

Description: This uses the csv's from the GEV_Analysis/ folder. The same calculation is performed on each
of the CF files and will be saved under a new name with the added column. For the purpose of this analysis CV is
calculated as the (upper - lower)/mean, the upper and lower bounds are the 90% CI.

"""

# Import Library
import gc
import os
import pandas as pd
from glob import glob

# Step 1: setup local environment and directory
# Prepare environment
pd.options.mode.copy_on_write = True
# get current working directroy
cwd_path = os.getcwd()
print(cwd_path)
# list files
cf_files = glob(cwd_path + f"/GEV_Analysis/*_CF_SSP*.csv", recursive=True)
print(cf_files)

# loop across the CF files to perform calculation
for i in range(len(cf_files)):
    print("starting with ", i+1, " : ", len(cf_files))

    # Step 2: load and read the file
    file = cf_files[i]
    cf_table = pd.read_csv(file)
    # create a copy of the data frame
    cf_table_copy = cf_table.copy(deep=True)

    # Step 3: Calculate the CVs of each station
    cf_table_copy['CV'] = (cf_table_copy['95th_percentile'] - cf_table_copy['5th_percentile']) / cf_table_copy['mean']

    # Step 4: Save as a new dataframe
    # create unique file saving name from original filename
    saving_file = file.replace('CF', 'CFcv')
    # save
    cf_table_copy.to_csv(saving_file, index=False)

    print("finished preparing", file, "and saved files")

    # Clean-up memory
    del cf_table, cf_table_copy
    gc.collect()

### FINISHED

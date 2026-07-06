# -*- coding: utf-8 -*-
"""
Date: 9/3/2025
Time: 4:27 PM
AIM: This script is writen to prepare GEV analysis precipitation events results (IDF estimate) for mapping interpolation

Description: GEV analysis results of 24-hr, 12-hr, and 6-hr precipitation events at stations are combined into separate
10-yr and 100-yr, and SSP2-4.5 and SSP5-8.5 dataframes while also ensuring there are no outliers remaining in the final
product for interpolation and maps.

2 Part code
    PART 1: Setup local variables and prepare for analysis
    PART 2: Prepare 24-hr, 12-hr, and 6-hr precipitation projected and observation IDFs

"""

import gc
import pandas as pd
import os
import datetime as dt
from glob import glob
import numpy as np
import pickle

##### PART 1: Setup local variables and prepare for analysis

def categorize_stations(x):
    """
    function to categorize the dataframe by station using grouping
    :param x: df.column
    :return : str
    """
    # create the list of grouped stations previously determined and used in the temporal disaggregation
    site_groups = dict(
        g1=['ELP', 'PEQ', 'GDP', 'FST', 'INK'],
        g2=['6R6', 'MAF', 'MRF', 'MDD', 'DUX', 'E38', 'ODO'],
        g3=['SNK', 'LRD', 'HRX', 'SJT', 'PYX', 'BPC', 'SWW', 'LBB', 'FTN', 'APY',
            'HBV', 'PVW', 'DHT', 'BPG', 'OZA', 'GNC', 'DRT', 'MFE', 'AMA'],
        g4=['UVA', 'BKS', 'ECU', 'JCT', 'COT', 'ALI', 'HRL', 'EBG', 'CDS', 'HHF', 'ABI'],
        g5=['BRO', 'COM', 'F05', 'PIL', 'T82', 'PEZ', 'AQO', 'BKD', 'SSF', 'NGP',
            'BEA', 'CWC', 'NQI', 'NOG', 'HDO', 'SPS', 'CVB', 'BMQ', 'DZB', 'RBO'],
        g6=['GRK', 'MKN', 'GOP', 'EDC', 'ERV', 'MNZ', 'RYW', 'LZZ', 'SEP', 'RPH',
            'RKP', 'HLR', 'GDJ', 'HYI', 'BAZ', 'CRP', 'MWL', 'SAT', 'ATT', 'AUS'],
        g7=['5C1', 'NFW', 'DFW', 'PWG', 'GYB', 'AFW', 'INJ', 'LUD', 'GPM',
            'RWV', 'FTW', 'JWY', 'DAL', 'DTO', 'XBP', 'GLE', 'FWS', 'GTU', 'RBD',
            'TPL', '0F2', 'LHB', 'ACT', 'GKY', 'RAS'],
        g8=['TRL', 'CLL', '11R', 'HQZ', '3T5', 'ADS', 'GYI', 'VCT', 'JDD', 'PKV', 'CRS', 'TKI', 'GVT'],
        g9=['GGG', 'SLR', 'PRX', 'JXI', 'TYR', 'PSX', 'PSN', 'JSO', 'GLS', 'OSA'],
        g10=['AXH', 'UTS', 'MCJ', 'SGR', 'IAH', 'OCH', 'RFI', 'DWH', 'CXO', 'LFK', 'BYY'],
        g11=['LVJ', 'HOU', '6R3', 'LBX', 'JAS'],
        g12=['ORG', 'BPT']
    )

    for category, keys in site_groups.items():
        if any(key in x for key in keys):
            return category
    return 'OTHER'


# Step 1: setup local environment and directory

# get current working directroy
cwd_path = os.getcwd()

# initial start time
ist = dt.datetime.now()
print("Initial script start time: ", ist)

# list files of interest
gev_files = glob(cwd_path + f"/GEV_Analysis/All_stations_*.csv", recursive=True)

# Step 2: load and prepare station descriptive data
# get ASOS data to add lat + lon + elev data
dir_master_ASOS = os.path.dirname(os.path.abspath(cwd_path)) + '/Data/ASOS/'
file_path_ASOS = os.path.join(dir_master_ASOS, 'TX_ASOS_Final.csv')
df_ASOS_stations = pd.read_csv(file_path_ASOS)
df_ASOS_stations = df_ASOS_stations[['stid', 'lat', 'lon', 'elev']]

# list of stations to remove because of their anomalous results
stids_removed = ['ARM', 'BGD', 'CPT', 'DLF', 'DYS', 'ILE', 'LNC', 'RND', 'SKF', 'T65']

# Step 3: setup dictionaries to retain details of what is done
# Dictionary for the outliers that were removed
dict_outlier_dfs = {}
# Dictionary for the data that remained
dict_dfs = {}

##### PART 2: Prepare 24-hr, 12-hr, and 6-hr precipitation projected and observation IDFs

# Prepare observation and bias corrected GCMs GEV results
for f in range(len(gev_files)):
    ist = dt.datetime.now()
    print("Initial script start time: ", ist)
    print("starting with ", f+1, " : ", len(gev_files))

    # Step 2: read GEV file and prepare base dataframe
    file = gev_files[f]
    orig_df = pd.read_csv(file, index_col=0)

    # add the lat, lon, and elevation station details
    df_join = pd.merge(orig_df, df_ASOS_stations, how="left", left_on="Station", right_on="stid",
                       validate="many_to_one")

    # Step 3: clean and correct GEV data
    # check to make sure to exclude the previously identified Stations that had clear error in data records
    df_clean = df_join.loc[~df_join['Station'].isin(stids_removed)].copy()

    # check for extra index column and if it exists remove it
    if 'Unnamed: 0' in df_clean.columns:
        print("yes")
        df_clean = df_clean.drop(columns=['Unnamed: 0'])

    # categorize the stations by similarity grouping
    df_clean['Groups'] = df_clean['Station'].apply(categorize_stations)
    # make sure that the necessary column names are all the same
    df_clean = df_clean.rename(columns={"Mean": "mean"})

    # OBSERVATIONAL DATA
    if "obs" in file:
        # Step 4: remove stations with unreasonable GEV results
        print(df_clean.groupby(['return_period'])['mean'].std())
        print("dataframe shape prior to outlier removal", df_clean.shape)
        # Calculate Z-score
        df_clean['z_score'] = df_clean.groupby(['return_period'], group_keys=False)['mean'].apply(
            lambda x: np.abs((x - x.mean()) / x.std()))
        df_NOoutliers = df_clean.loc[lambda df: df['z_score'] <= 3].copy()
        df_outliers = df_clean.loc[lambda df: df['z_score'] > 3].copy()

        # add the outliers removed to the dictionary
        dict_outlier_dfs[file[file.rfind('All'): -4]] = df_outliers
        dict_dfs[file[file.rfind('All'): -4]] = df_clean
        print("dataframe shape after outlier removal", df_NOoutliers.shape)

        # Step 5: create new 10-yr and 100-yr dataframes
        # filter and pull data by scenario and return period
        df_10 = df_NOoutliers.loc[(df_NOoutliers['return_period'] == 0.9)].copy().drop(columns= ['z_score', 'stid'])
        df_100 = df_NOoutliers.loc[(df_NOoutliers['return_period'] == 0.99)].copy().drop(columns= ['z_score', 'stid'])
        print("all dataframes created")

        # Step 6: Save these dataframes
        # create unique file saving name from original filename
        saving_file1 = file.replace('summary', '10yr')
        saving_file2 = file.replace('summary', '100yr')
        # save
        df_10.to_csv(saving_file1, index=False)
        df_100.to_csv(saving_file2, index=False)
        print("finished preparing", file, "and saved files")
    # BIAS CORRECTED GCM DATA
    else:
        # Step 7: remove stations with unreasonable GEV results
        print(df_clean.groupby(['Scenario', 'return_period'])['mean'].std())
        print("dataframe shape prior to outlier removal", df_clean.shape)

        # Calculate Z-scores of mean
        df_clean['z_score'] = df_clean.groupby(['Scenario', 'return_period'], group_keys=False)['mean'].apply(
            lambda x: np.abs((x - x.mean()) / x.std()))

        # remove outliers
        df_NOoutliers = df_clean.loc[lambda df: df['z_score'] <= 3].copy()
        # pull outliers
        df_outliers = df_clean.loc[lambda df: df['z_score'] > 3].copy()

        # add the outliers removed to the dictionary
        dict_outlier_dfs[file[file.rfind('All'): -4]] = df_outliers
        dict_dfs[file[file.rfind('All'): -4]] = df_clean
        print("dataframe shape after outlier removal", df_NOoutliers.shape)

        # Step 8: create new 10-yr and 100-yr dataframes
        # filter and pull data by scenario and return period
        df_ssp245_10 = df_NOoutliers.loc[(df_clean['Scenario'] == "SSP245") &
                                         (df_NOoutliers['return_period'] == 0.9)].copy().drop(
            columns= ['z_score', 'stid'])
        df_ssp585_10 = df_NOoutliers.loc[(df_clean['Scenario'] == "SSP585") &
                                         (df_NOoutliers['return_period'] == 0.9)].copy().drop(
            columns= ['z_score', 'stid'])
        df_ssp245_100 = df_NOoutliers.loc[(df_clean['Scenario'] == "SSP245") &
                                          (df_NOoutliers['return_period'] == 0.99)].copy().drop(
            columns= ['z_score', 'stid'])
        df_ssp585_100 = df_NOoutliers.loc[(df_clean['Scenario'] == "SSP585") &
                                          (df_NOoutliers['return_period'] == 0.99)].copy().drop(
            columns= ['z_score', 'stid'])
        print("all dataframes created")

        # Step 9: Save these dataframes
        # create unique file saving name from original filename
        saving_file1 = file.replace('summary', 'SSP245_10yr')
        saving_file2 = file.replace('summary', 'SSP585_10yr')
        saving_file3 = file.replace('summary', 'SSP245_100yr')
        saving_file4 = file.replace('summary', 'SSP585_100yr')
        # save
        df_ssp245_10.to_csv(saving_file1, index=False)
        df_ssp585_10.to_csv(saving_file2, index=False)
        df_ssp245_100.to_csv(saving_file3, index=False)
        df_ssp585_100.to_csv(saving_file4, index=False)
        print("finished preparing", file, "and saved files")

        # clean-up memory
        del df_ssp245_100, df_ssp245_10, df_ssp585_10, df_ssp585_100, (
            saving_file1), saving_file2, saving_file3, saving_file4
    t2 = dt.datetime.now()
    print("loop time length: ", t2 - ist)
    # clean-up memory
    del file, orig_df, df_join, df_clean,df_NOoutliers, df_outliers
    gc.collect()

# inspect the station outliers that were removed
print("Station outlier data removed from their datasets")
print(dict_outlier_dfs)

# save the dictionary of outliers removed
with open("dict_outliers_removed.pkl", "wb") as file:
    pickle.dump(dict_outlier_dfs, file)

### FINISHED
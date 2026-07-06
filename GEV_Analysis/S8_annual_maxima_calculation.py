"""
Author: samantha murray
Date: 9/3/2025
Time: 1:44 PM
AIM: This script is writen to compute the extreme precipitation as annual maxima for observation and bias corrected GCMs

Description: Calculate the annual maximum to compute the extreme precipitation distribution analysis for 24-HR, 12-HR,
 and 6-HR precipitation. Script built on the framework of a multi-model analysis.
 2 Part code to calculate extreme precipitation values based on data structures
    PART 1: calculate 24-HR annual maximum bias corrected projected GCMs and observation precipitation
    PART 2: calculate 12-HR & 6-HR annual maximum bias corrected projected GCMs precipitation

"""
import gc
import xarray as xr
import numpy as np
import pandas as pd
import os
import math
import datetime as dt
import glob
from xclim import ensembles
from xclim import units
from functools import reduce

##### PART 1: calculate 24-HR annual maximum bias corrected projected GCMs and observation precipitation

### Function to get the bias corrected data filepaths from current directory
# input parameters: directory pathway to search through; station_id the ASOS station ID of interest
# Example usage
# directory = cwd_path
# station_id = 'ADS'
def list_bc_csvs(directory, station_id):
    """
    Function to get and list the bias corrected data filepaths from directory pre-developed pathway and cleans away
    excess files
    input parameters:
    directory == pathway to search through,
    station_id == the ASOS station ID of interest
    """
    all_files = glob.glob(directory + f"BiasCorrection/BiasCorrected_Data/*_biascorr/**/{station_id}_*.csv",
                          recursive=True)
    return all_files


# function to create saving folder pathway for each model and scenario if it does not exist
# input parameter: folder pathway you want to save to
def create_saving_folders(folder_path):
    """
    Function to create saving folder pathway for each if it does not exist
    excess files
    input parameters:
    folder_path == directory to create and its pathway
    """
    if not os.path.exists(folder_path):
        # if here current thread is stopped and the same dir is created in other thread
        # the next line will raise an exception
        os.makedirs(folder_path)
        print("directory made")
    else: print("directory already exists")

# Function to load Observation data file path by station ID and prepare to match data structure of bias corrected dataframes
# input parameters: ASOS netcdf file path, ASOS dictionary, station id, variable name extracting from .nc
# Example usage
# station_id = 'DAL'
# variable_OBS = 'pr' (i.e., precipitation)
def get_obs_station_data(directory, station_id):
    """
    Function to load Observation data file path by station ID and prepare to match data structure of
    bias corrected dataframes. Pulls and uses the cleaned Observation data.
    input parameters:
    directory == current working root directory 
    station_id == station ID; string
    """
    # locate the dictionary start date for station
    daily_file = (directory + f"Data/ASOS/ASOS_Stations_Data_Preparation/PreparedData/Daily/"
                              f"{station_id}_daily_precipitation_cleaned.csv")
    # read file
    df = pd.read_csv(daily_file)
    # rename precipitation column
    df = df.rename(columns={'Precipitation': 'pr'})
    # set column as date time objects
    df['time'] = pd.to_datetime(df['time'])
    # Keep only the observations after 1950
    df = df[df['time'] >= '1950-01-01']
    # create new date columns: day of year, month, year
    df['doy'] = df['time'].dt.dayofyear
    df['month'] = df['time'].dt.month
    df['year'] = df['time'].dt.year
    # set index
    df = df.set_index('time')
    # now restructure this df to resemble the structure of observation df
    # add columns to match
    df['type'] = 'observation'
    df['scenario'] = 'observation'
    df['model'] = 'observation'
    df['method'] = 'observation'
    return df


# function to read listed filepaths from current directory and adds all dataframes into one dataframe
# use with the files gathered from the list_bc_csvs function
def process_biascorr_df(files):
    """
    Function to load and prepare Bias Corrected GCM files. Create columns related to date, station ID, model,
    scenario, method, and type for complete data description.
    Used with the files gathered from the list_bc_csvs function
    input parameters:
    files == bias corrected data filepaths from current directory at a ASOS Station
    """
    all_data = [] # empty list to append dataframes to
    for filename in files:
        # read file
        df = pd.read_csv(filename)
        # set datetime
        df['time'] = pd.to_datetime(df['time'])
        # create new date columns: day of year, month, year
        df['doy'] = df['time'].dt.dayofyear
        df['month'] = df['time'].dt.month
        df['year'] = df['time'].dt.year
        # set the index
        df = df.set_index('time')
        # get string index to use for unique identifiers
        str_i = filename.find(
            '5\\') + 2  # '5\\' is unique string found at the end of all files before naming convention
        ## name string to use
        new_name = filename[str_i:-4]
        # correct names to remove information that won't be needed for classifications
        new_name = new_name.replace('_2015_', '_')
        new_name = new_name.replace('_nq50', '_')
        ## rename
        df['identifier'] = new_name
        # split the identifier column by string
        # *Note all stings follow the format "stid_model_scenario_methodtype_methodP1_P2_P3"
        # example:"RBD_GFDL-ESM4_SSP585_nq50multiplication_quantile_mapping"
        # df[['stid', 'model', 'scenario', 'type', 'methodP1', 'methodP2', 'methodP3']] = (
        #     df['identifier'].str.split('_', expand=True))
        name_split = new_name.split('_')
        df['stid'] = name_split[0]
        df['model'] = name_split[1]
        df['scenario'] = name_split[2]
        df['type'] = name_split[3]
        df['methodP1'] = name_split[4]
        df['methodP2'] = name_split[5]
        if len(name_split) > 6:
            df['methodP3'] = name_split[6]
            df['method'] = df['methodP1'] + ' ' + df['methodP2'] + ' ' + df['methodP3']
        else:
            # now concat all methodP columns into one 'method' column
            df['method'] = df['methodP1'] + ' ' + df['methodP2']
        df['method'] = df['method'].str.strip()
        # drop the old columns
        df = df.drop(columns=df.columns[9:-1])
        df = df.drop(columns='identifier')
        # append
        all_data.append(df)
    final_df = pd.concat(all_data, axis=0, ignore_index=False)
    return final_df


# Function to calculate annual 24-HR maximums dataframes for each ASOS station
def create_24hr_maximums(
        df: pd.DataFrame,
        is_bc_df: str):
    # Create a copy of the original DataFrame to avoid modifying it
    new_df = df.copy()
    if is_bc_df in ['y', 'ye', 'yes']:
        # reduce the methods section to include methods and type
        new_df["methods"] = new_df["type"] + ' ' + new_df["method"]
        # drop old columns
        new_df = new_df.drop(columns=['type', 'method'])
        # calculate the daily max for each scenario, model and method
        df_annual_24hrmax = new_df.groupby(['stid', 'model', 'scenario', 'methods', 'year'])[['pr']].max()
        print('completed 24-HR maximums of bias corrected data')
        return df_annual_24hrmax
    if is_bc_df in ['n', 'no', 'nop', 'nope']:
        # calculate the daily max for each scenario, model and method
        df_annual_24hrmax = new_df.groupby(['stid', 'scenario', 'year'], as_index=False)[['pr']].max()
        print('completed 24-HR maximums of observational data')
        return df_annual_24hrmax


# initial start time
ist = dt.datetime.now()
print("Initial script Part 1 start time: ", ist)

# get current working directory
cwd_path = os.getcwd()

# Step 1: load and prepare ASOS station data
# load ASOS data
asos_xy = pd.read_csv(cwd_path + '/Data/ASOS/TX_ASOS_Final.csv')
# Adjust display settings to show all columns
pd.set_option('display.max_columns', 10)
# get station id's found in asos_final
sites = asos_xy['stid'].values.tolist()
# set stid as index
asos_xy = asos_xy.set_index('stid')

# Step 2: set loop parameters to loop through GCMs and scenarios
# climate scenarios
SSPs = ['245', '585']
# models
models = ['CanESM5', 'GFDL-ESM4', 'BCC-CSM2-MR', 'EC-Earth3', 'CNRM-CM6-1', 'ACCESS-CM2', 'CNRM-ESM2-1']

# create a list for a new dataframe to hold the determined Threshold value
dfs = []

# create saving folders
create_saving_folders(cwd_path + "Station_combined_data")
create_saving_folders(cwd_path + 'GEV_Analysis/Station_24HR_Maximums')
create_saving_folders(cwd_path + 'GEV_Analysis/Station_24HR_Maximums/Observation_24HRmax')
create_saving_folders(cwd_path + 'GEV_Analysis/Station_24HR_Maximums/Station_biascorr_24HRmax')
create_saving_folders(cwd_path + 'GEV_Analysis/Station_12HR_Maximums/Station_biascorr_12HRmax')
create_saving_folders(cwd_path + 'GEV_Analysis/Station_6HR_Maximums/Station_biascorr_6HRmax')

# Prepare bias corrected data into 1 csv and get 24 hr max
# loop across station IDs
for si in range(len(sites)):
    print("starting with ", si, " : ", len(sites))
    # define station id
    stid = sites[si]
    print(stid)
    # variable of interest
    var = 'pr'

    # Step 3: get station bias corrected climate csvs
    # get bias corrected files for station
    bias_corr_files = list_bc_csvs(cwd_path, stid)
    # check length
    print("# of files:", len(bias_corr_files))

    # bypass if there are no bias corrected files for this station
    if len(bias_corr_files) == 0:
        pass
    else:
        # Step 4: Load and Prepare bias corrected dfs
        stid_bias_corr_df = process_biascorr_df(bias_corr_files)

        ## Check for negative values of the bias correction results
        negative_pr = stid_bias_corr_df[stid_bias_corr_df['pr'] < 0]
        print(negative_pr.describe(include=object))
        print('Determine which methods have negative values:', negative_pr['method'].unique())
        print('Determine which types have negative values:', negative_pr['type'].unique())
        # it is expected to have negative values for additive effects of linear scaling method

        # Step 5: Load and Prepare station observation data
        stid_obs_df = get_obs_station_data(directory=cwd_path, station_id=stid)

        # reorganize columns
        cols_org = list(stid_bias_corr_df.columns)
        stid_obs_df = stid_obs_df[cols_org]

        # Step 6: prepare to merge all these dataframes into one large one for seaborn plotting
        print("Check that these two dataframes will match for future uses")
        # double check all these dataframes are compatible
        print("set difference:", np.setdiff1d(stid_bias_corr_df.columns, stid_obs_df.columns))
        print("Column Names of Bias Corrected data: ", stid_bias_corr_df.columns)
        print("Column Names of Observation_results data: ", stid_obs_df.columns)
        col_names = list(stid_bias_corr_df.columns)
        # now we know they match we can continue

        # Save
        stid_bias_corr_df.to_csv(f'Station_combined_data/{stid}_bias_corr_df.csv', index=True)
        stid_obs_df.to_csv(f'Station_combined_data/{stid}_observation_df.csv', index=True)

        # Step 8: Generate annual 24-HR maximums for Bias Corrected data

        # bias corrected df
        df_annual_24hrmax_bc = create_24hr_maximums(df=stid_bias_corr_df, is_bc_df='y')
        # reset index
        df_annual_24hrmax_bc = df_annual_24hrmax_bc.reset_index()
        # Save
        df_annual_24hrmax_bc.to_csv(
            f'GEV_Analysis/Station_24HR_Maximums/Station_biascorr_24HRmax/{stid}_bias_corr_24HRmax_df.csv',
            index=False)

        # Step 9: Generate annual 24-HR maximums for Observational data

        # observational df
        df_annual_24hrmax_obs = create_24hr_maximums(df=stid_obs_df, is_bc_df='n')
        # reset index
        df_annual_24hrmax_obs = df_annual_24hrmax_obs.reset_index()
        # Save
        df_annual_24hrmax_obs.to_csv(
            f'GEV_Analysis/Station_24HR_Maximums/Observation_24HRmax/{stid}_obs_24HRmax_df.csv',
            index=False)

        # clean environment
        del stid_obs_df, stid_bias_corr_df, bias_corr_files, negative_pr, df_annual_24hrmax_obs, df_annual_24hrmax_bc
        gc.collect()

# ending time
et = dt.datetime.now()
print("Script Part 1 end time: ", et)
print("Part 1 execution time: ", et - ist)

##### PART 2: calculate 12-HR & 6-HR annual maximum bias corrected projected GCMs precipitation

# Function to read listed filepaths from current directory and adds all dataframes into one dataframe
# input parameters: dask dataframe, list of csv files
# use with the files gathered from the list_bc_csvs function
def process_biascorr_pred_df(files, station_id):
    all_data = [] # empty list to append dataframes to
    for filename in files:
        # read file
        df = pd.read_csv(filename, index_col=0)
        # remove the un-needed columns (i.e., model inputs)
        df = df.drop(columns=['elev', 'pr_monthly_T', 'pr_annual_T', 'pr_1day_T', 'pr_+1day_T', 'pr_-1day_T',
                              'pr_+1monthly_T', 'pr_-1monthly_T', 'Diff_6hr', 'Diff_12hr'])
        # set datetime
        df['Date'] = pd.to_datetime(df['Date'])
        # create new date columns: month, year
        df['year'] = df['Date'].dt.year
        # set the index
        df = df.set_index('Date')
        # get string index to use for unique identifiers
        str_i = filename.find(
            station_id) + 4  # '5\\' is unique string found at the end of all files before naming convention
        ## name string to use
        new_name = filename[str_i:-4]
        # correct names to remove information that won't be needed for classifications
        new_name = new_name.replace('_2015_', '_')
        new_name = new_name.replace('_nq50', '_')
        ## rename
        df['identifier'] = new_name
        # split the identifier column by string
        # *Note all stings follow the format "stid_model_scenario_methodtype_methodP1_P2_P3"
        # example:"GFDL-ESM4_SSP585_nq50multiplication_quantile_mapping"
        # df[['model', 'scenario', 'type', 'methodP1', 'methodP2', 'methodP3']] = (
        #     df['identifier'].str.split('_', expand=True))
        name_split = new_name.split('_')
        df['model'] = name_split[0]
        df['scenario'] = name_split[1]
        df['type'] = name_split[2]
        df['methodP1'] = name_split[3]
        df['methodP2'] = name_split[4]
        if len(name_split) > 5:
            df['methodP3'] = name_split[5]
            df['method'] = df['methodP1'] + ' ' + df['methodP2'] + ' ' + df['methodP3']
        else:
            # now concat all methodP columns into one 'method' column
            df['method'] = df['methodP1'] + ' ' + df['methodP2']
        df['method'] = df['method'].str.strip()
        # drop the old columns
        df = df.drop(columns=df.columns[11:-1])
        df = df.drop(columns='identifier')
        # append
        all_data.append(df)
    final_df = pd.concat(all_data, axis=0, ignore_index=False)
    return final_df

# Function to calculate annual 12- & 6-HR maximums dataframes for each ASOS station
def create_annual_maximums_bc(
        # station_id: str,
        df: pd.DataFrame):
    # Create a copy of the original DataFrame to avoid modifying it
    new_df = df.copy()
    # reduce the methods section to include methods and type
    new_df["methods"] = new_df["type"] + ' ' + new_df["method"]
    # drop old columns
    new_df = new_df.drop(columns=['type', 'method'])
    # calculate the 6-HR max for each scenario, model and method
    df_annual_6hrmax = new_df.groupby(['stid', 'model', 'scenario', 'methods', 'year'])[['Predicted_6hr_pr']].max()
    # reset index
    df_annual_6hrmax = df_annual_6hrmax.reset_index()
    # # Save
    # df_annual_6hrmax.to_csv(
    #     f'Station_6HR_Maximums/Station_biascorr_6HRmax/{station_id}_bias_corr_6HRmax_df.csv', index=False)
    print('completed 6-HR maximums of bias corrected data')
    # calculate the 12-HR max for each scenario, model and method
    df_annual_12hrmax = new_df.groupby(['stid', 'model', 'scenario', 'methods', 'year'])[['Predicted_12hr_pr']].max()
    # reset index
    df_annual_12hrmax = df_annual_12hrmax.reset_index()
    # # Save
    # df_annual_12hrmax.to_csv(
    #     f'Station_12HR_Maximums/Station_biascorr_12HRmax/{station_id}_bias_corr_12HRmax_df.csv', index=False)
    print('completed 12-HR maximums of bias corrected data')
    return df_annual_6hrmax, df_annual_12hrmax


# initial start time
ist2 = dt.datetime.now()
print("Initial script Part 2 start time: ", ist2)

# get current working directory
cwd_path = os.getcwd()

# Step 1: load and prepare ASOS station data
# load ASOS data
asos_xy = pd.read_csv(cwd_path + '/Data/ASOS/TX_ASOS_Final.csv')
# get station id's found in asos_final
sites = asos_xy['stid'].values.tolist()
# set stid as index
asos_xy = asos_xy.set_index('stid')

# Step 2: set loop parameters to loop through station IDs

# create a list for a new dataframe to hold the determined Threshold value
dfs = []

# Prepare bias corrected temporally disaggregated events to get the annual maximum for 12-hr and 6-hr events
for si in range(len(sites)):
    print("starting with ", si+1, " : ", len(sites))
    # Step 3: get station bias corrected climate csvs
    
    # define station id
    stid = sites[si]
    print(stid)
    
    # get bias corrected files for station
    bias_corr_pred_files = glob.glob(cwd_path + f"/Temporal disaggregation/model_results/**{stid}_*.csv",
                                     recursive=True)
    # check length
    print("# of files:", len(bias_corr_pred_files))
    # bypass if there are no bias corrected files for this station
    if len(bias_corr_pred_files) == 0:
        pass
    else:
        # Step 4: Load and Prepare bias corrected dfs
        pred_stid_bias_corr_df = process_biascorr_pred_df(bias_corr_pred_files, station_id=stid)
        # Save
        pred_stid_bias_corr_df.to_csv(f'Station_combined_data/{stid}_bias_corr_6hr-12hr_df.csv', index=True)

        # Step 5: Calculate annual maximums of bias corrected GCMs 12-HR and 6-HR precipitation
        # 6-HR, 12-HR
        df_amax_6hr, df_amax_12hr = create_annual_maximums_bc(df= pred_stid_bias_corr_df)
        # Save
        df_amax_6hr.to_csv(
            f'GEV_Analysis/Station_6HR_Maximums/Station_biascorr_6HRmax/{stid}_bias_corr_6HRmax_df.csv',
            index=False)
        print("6-HR annual maximum table saved")
        df_amax_12hr.to_csv(
            f'GEV_Analysis/Station_12HR_Maximums/Station_biascorr_12HRmax/{stid}_bias_corr_12HRmax_df.csv',
            index=False)
        print("12-HR annual maximum table saved")

        # clean environment
        del df_amax_6hr, df_amax_12hr, stid, pred_stid_bias_corr_df, bias_corr_pred_files
        gc.collect()


# ending time
et = dt.datetime.now()
print("Script Part 2 end time: ", et)
print("Part 2 execution time: ", et - ist2)
print("FINISHED")

### FINISHED
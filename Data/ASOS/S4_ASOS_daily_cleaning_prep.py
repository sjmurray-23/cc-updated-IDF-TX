# -*- coding: utf-8 -*-

""""
AIM: Script for preparing precipitation to apply bias correction techniques and temporal disaggregation
Cleaning and preparing previously extracted and generated ASOS station data for analysis
Data from Part 4 daily precipitation generation of ASOS_data_generation.py is cleaned
Data from Part


# 4 Part code for preparing previously extracted ASOS station data for analysis

# PART 1: clean missing data in ASOS station daily data
# PART 2: identify and exclude ASOS stations with anomalies
# PART 3: Prepare Observation daily for netcdf file
# PART 4: generate Observation daily netcdf file

"""""

import pandas as pd
from glob import glob
from os import path
from os import getcwd
import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
import datetime
import xarray as xr
import gc
import datetime as dt

# setup
cwd = getcwd()
data_folder = r'Data\ASOS\PreparedData\Daily'
daily_files = glob(cwd + data_folder + "/*.csv")

# PART 1: clean ASOS station daily data
# Aim: Locate and remove missing data for each ASOS station

# Create a dictionary to store the different details of each station
ASOS_cleaning_dic = {}
# Ending dictionary structure
# ASOS_cleaning_dic = {
#     'stid': 'string', # Station ID
#     'starting_observations': 'integer', # number of observations before cleaning
#     'ending_observations' : 'integer', # number of observations after cleaning
#     'Years' : 'list', # all the years of data before cleaning
#     'Yearly_precip' : 'list', # yearly precipitation totals before cleaning
#     'Corrected' : 'string' # Y/N if the data was corrected for large (>365 days) data gaps
#     'new_start_date': 'string', # new starting date after cleaning
# }


# Prepare cleaned daily precipitation dataframes by looping through files
for fi in range(len(daily_files)):
    print(fi, "out of", len(daily_files))

    # Step 1: read file
    f = daily_files[fi]
    daily = pd.read_csv(f)
    print(daily.head())

    # Step 2: get stid name and total number of observations
    site = daily['stid'].unique()[0]
    print(site)
    s_obs = len(daily.index)

    # Step 3: format time
    daily['time'] = pd.to_datetime(daily['time'])
    print(daily.info())
    # set time as index
    daily.set_index('time', inplace=True)

    # Step 4: identify the gaps in data (numerous consecutive 0 values of annual rainfall)
    # generate yearly rainfall
    yearly = daily.resample('YE').sum()
    print(yearly.head())
    # correct the stid column
    yearly['stid'] = site
    years = yearly.index.year
    preci = yearly['Precipitation'].values
    bool_preci = (preci == 0)

    # Step 5: remove incomplete records of data before the date of data lost so the data set is continuous
    if sum(bool_preci) > 1:
        print(f'{site} did have to be cleaned')
        # dictionary key of cleaned
        c = 'Y'
        ## identify the last year when there is no data
        zero_yr = yearly[yearly.Precipitation == 0]
        last_yr = zero_yr.index.max()
        zero_daily = daily.loc[:last_yr]
        # get the difference between two dataframes
        daily_cleaned = daily[~daily.index.isin(zero_daily.index)]
    else:  # simply rename the file
        # dictionary key of cleaned
        c = 'N'
        daily_cleaned = daily
    print('Completed cleaning process')

    # Step 6: fill the dictionary
    # ending number of observations
    e_obs = len(daily_cleaned.index)
    # get new start date
    start = daily_cleaned.index.min()
    temp_dic = {
        site: {
            'stid': site,
            'starting_observations': s_obs,
            'ending_observations': e_obs,
            'Years': years,
            'Yearly_precip': preci,
            'Corrected': c,
            'new_start_date': start}
    }
    ASOS_cleaning_dic.update(temp_dic)
    print('added to dictionary')
    del temp_dic

    # Step 7: plot
    daily_cleaned.plot(kind='line', color='blue', linestyle='-', linewidth=2)
    plt.xlabel('Date')
    title = site + ' Daily Precipitation - Cleaned'
    plt.ylabel('Daily Precipitation (in)')
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("fig/" + title + ".png")

    # Step 8: save the file
    outfile = site + "_daily_precipitation_cleaned.csv"
    daily_cleaned.to_csv(path.join(r"PreparedData\Daily", outfile))
    print(f'{site} is complete and saved')

# Save dictionary to a pickle file
with open(data_folder + 'ASOS_daily_cleaning_dict.pkl', 'wb') as fp:
    pickle.dump(ASOS_cleaning_dic, fp)
    print('dictionary saved successfully to file')

# PART 2: identify and exclude ASOS stations with anomalies
# Aim: review data dictionary created and identify the stations that:
# do not have enough measures or at least 6 years of data after being cleaned
# too little rainfall or yearly rainfall <1 or no precipitation values (0)

# Setup
cwd = getcwd()
print(cwd)
data_folder = r'Data\ASOS\PreparedData\Daily'

# Step 1: load data

ASOS_stations = pd.read_csv('Data\ASOS\TX_ASOS_Stations_CleanClipped.csv')
# get the station IDs, lat, lons
ASOS_stations = ASOS_stations[['stid', 'lat', 'lon']]
# get station id's found in asos_final
sites = ASOS_stations['stid'].values.tolist()
# set stid as index
ASOS_stations.set_index('stid', inplace=True)

# open dictionary
with open(data_folder+'ASOS_daily_cleaning_dict.pkl', 'rb') as file:
    ASOS_dict = pickle.load(file)

# Step 2: evaluate the dictionary
# check the number of stations uploaded to the dictionary
keys = list(ASOS_dict.keys())  # save as a list to iterate through later
print(len(keys))
# check that the keys will access the dictionary and the structure of it
print(ASOS_dict[keys[0]])

# identify which stations were corrected
# empty list to save to
stids_were_corrected = []
for key in keys:
    corr = ASOS_dict[key].get('Corrected')
    if corr == 'Y':
        stids_were_corrected.append(ASOS_dict[key].get('stid'))

# the number of the stations that were corrected
len(stids_were_corrected)

# Step 3: Identify any stations that had only odd precipitation values (<1) or no precipitation values (0)
# create empty list to save to
stids_with_odd_precip = []
for key in keys:
    precip = ASOS_dict[key].get('Yearly_precip')
    # using <1 for all the summed yearly precipitation values, since all years with no/small precipitation cannot be analyzed or trusted as accurate
    if all(precip < 1):
        stids_with_odd_precip.append(ASOS_dict[key].get('stid'))

# the number of stations with anamolous precipitation, figure charts show that 19 stations have no values
len(stids_with_odd_precip)
print(stids_with_odd_precip)

# 19 stations that figures show no values
stids_figs_check = ['TME', 'T69', 'SPL', 'SOA', 'SHP', 'RPE', 'NMT', 'MIU', 'LBR', 'HQI', 'HHV', 'F12', 'EMK', 'EFD',
                    'CNW', 'BQX', 'BBF', 'BMT', '25T']

# check these keys
for key in stids_figs_check:
    print(ASOS_dict[key])
# many of the above stations ended with 0 observations

# Step 4: Identify any stations ended with 0 observations
# empty list to save to
stids_no_ending = []
for key in keys:
    end = ASOS_dict[key].get('ending_observations')
    if end == 0:
        stids_no_ending.append(ASOS_dict[key].get('stid'))
print(stids_no_ending)
# SHP was not in this list, 18 in total
# checking back there is only 1 year of observations that = 0

# Step 5: Identify stations have at least 6 years of data after being cleaned
# create empty list to save to
stids_less_than_6y = []
for key in keys:
    y = len(ASOS_dict[key].get('Years'))
    if y < 6:
        stids_less_than_6y.append(ASOS_dict[key].get('stid'))
len(stids_less_than_6y)
# inspect these
# check these keys
for key in stids_less_than_6y:
    print(ASOS_dict[key])
    print(len(ASOS_dict[key].get('Years')))

# Step 6: Combine these stations to create a file of stations to remove and cleaned asos id file
# identify if there is a difference between dictionary keys list and the ASOS station clean clipped list
difference1 = list(set(keys) - set(sites))
print(difference1)
print(len(difference1))
# remove these from the identified stations to remove lists to resolve errors
# stations with odd precipitation values
stids_with_odd_precip_cln = list(set(stids_with_odd_precip) - set(difference1))
print(stids_with_odd_precip_cln)
print(len(stids_with_odd_precip_cln))
# stations with no ending rainfall
stids_no_ending_cln = list(set(stids_no_ending) - set(difference1))
print(stids_no_ending_cln)
print(len(stids_no_ending_cln))
# stations with less than 6 years of data
stids_less_than_6y_cln = list(set(stids_less_than_6y) - set(difference1))
print(stids_less_than_6y_cln)
print(len(stids_less_than_6y_cln))

# combine these lists and remove the duplicates to get one master list of station IDs that will be removed from the analysis
stids_to_drop = stids_with_odd_precip_cln + stids_no_ending_cln + stids_less_than_6y_cln
print(len(stids_to_drop))
# convert sets to one list also removing duplicates
stids_to_drop = list(set(stids_to_drop))
print('list after removing duplicates ', len(stids_to_drop))
print(stids_to_drop)

# Step 8: Identify stations to remove stations that are temporally too short that did not remove properly before
# Filter stations based on start date
start_date_too_early = [stid for stid in sites if ASOS_dict[stid].get('new_start_date') > '2010-12-31']
# inspect
print(len(start_date_too_early))

# Step 9: Identify any stations that had too much precipitation in a year (>100)
# create empty list to save to
stids_too_much_precip = []
# use dictionary to find these stations
for key in keys:
    precip = ASOS_dict[key].get('Yearly_precip')
    # using <1 for all the summed yearly precipitation values, since all years with no/small precipitation cannot be analyzed or trusted as accurate
    if any(precip > 100):
        stids_too_much_precip.append(ASOS_dict[key].get('stid'))

# the number of stations with anamolous precipitation, figure charts show that 19 stations have no values
print(len(stids_too_much_precip))
print(stids_too_much_precip)


# Step 10: Remove stations from the ASOS dataframe
# join the lists with the other stations to remove
stids_to_drop = stids_to_drop + start_date_too_early + stids_too_much_precip

# save this list
with open('ASOS_stationIDs_to_remove.txt', 'w') as file:
    for item in stids_to_drop:
        file.write(f"{item}\n")

# remove list of stations
remaining_sites = list(set(sites) - set(stids_to_drop))
print(len(remaining_sites))

ASOS_stations_clean = ASOS_stations.loc[remaining_sites]

# append the new start date to the ASOS station
ASOS_stations_clean['new_start_date'] = [ASOS_dict[stid].get('new_start_date') for stid in ASOS_stations_clean.index]
# double check
print(ASOS_stations_clean)
print(ASOS_dict['GGG'].get('new_start_date'))
print(ASOS_dict['NFW'].get('new_start_date'))

# Save the new dataframe
ASOS_stations_clean = ASOS_stations_clean.reset_index()
ASOS_stations_clean.to_csv(r'Data\ASOS\TX_ASOS_locs_Final.csv', index=False)

# PART 3: Prepare ASOS daily for netcdf file
# Aim: prepare to match units of GCMs and combine ASOS observation daily data into a single dataframe
# to be converted to a netcdf

# Step 1: load csv file
# Load ASOS station locations
# Final has only the stations that will be used in the analysis i.e., clean data
asos_loc_df = pd.read_csv(r'Data\ASOS\TX_ASOS_locs_Final.csv')
# Load cleaned daily data
data_folder = 'Data/ASOS/PreparedData/Daily'
daily_files = glob(data_folder + "/*cleaned.csv")

# Step 2: subset station files to only final stations identified
# final station id's found in asos_final
sites = asos_loc_df['stid'].values.tolist()
daily_files_clean = [file for file in daily_files if any(s in file for s in sites)]

# Step 3: Prepare station data and join to single list
# create empty list
asos_dfs = []
# loop to read and append dataframes
for fi in range(len(daily_files_clean)):
    print(fi, "out of", len(daily_files_clean))
    f = daily_files_clean[fi]

    # read csv
    preci = pd.read_csv(f)
    print(preci.head())

    # convert to datetime object
    preci['time'] = pd.to_datetime(preci['time'])
    # set time as index
    preci.set_index('time', inplace=True)
    # make sure the dataframe is sorted by datetime index
    preci.sort_index()

    # convert precipitation from inches to mm to match the GCM units
    preci['Precipitation'] = preci['Precipitation'] * 25.4

    # reshape data
    site = preci['stid'].unique()[0]
    preci = preci.rename(columns={'Precipitation': site})
    preci = preci.drop(columns=['stid'])

    # append to list
    asos_dfs.append(preci)
    print("finished preparing and joining", site)

# Step 4: Combine these into one big table
preci_daily_df = pd.concat(asos_dfs, axis=1)
# save this dataframe
outfile = "All_ASOS_daily_precipitation.csv"
preci_daily_df.to_csv(path.join(data_folder, outfile))

# Step 5: further prepare dataframe time and date to match the LOCA2 GCMs
# convert 'time' to datetime object
preci_daily_df['time'] = pd.to_datetime(preci_daily_df['time'])
# correct datetime to match Climate models by adding offset
offset = pd.tseries.frequencies.to_offset("12h")
preci_daily_df['time'] = preci_daily_df['time'] + offset

# Step 6: reshape dataframe to make station id column ["stid"] and precipitation values into another ["pr"]
# pull station ids
sites = preci_daily_df.columns[1:]
# reshape
preci_daily_dfmelt = preci_daily_df.melt(id_vars=['time'], var_name="stid", value_name="pr")
# set the indexes
preci_daily_dfind = preci_daily_dfmelt.set_index(['time', 'stid'])

# convert precip to kg m-2 s-1
preci_daily_dfind["pr"] = preci_daily_dfind["pr"] / 86400

# Step 7: Apply final corrections
# if applicable, drop the additional stations found to have spatial overlap and to be the less complete dataset
asos_loc_df = asos_loc_df[asos_loc_df.stid != 'T89']
# save this dataframe
outfile = r'Data\ASOS\TX_ASOS_locs_Final.csv'
asos_loc_df.to_csv(outfile)

# get the lat lons to combine with precip data
asos_xy_df = asos_loc_df[['stid', 'lat', 'lon']]
# merge the two tables to give them
asos_pr_daily_loc = pd.merge(asos_xy_df, preci_daily_dfind.reset_index(), on="stid", how="inner",
                             validate="one_to_many")

# save this dataframe
outfile = "All_ASOS_daily_precipitation_Long.csv"
asos_pr_daily_loc.to_csv(path.join(data_folder, outfile))

# PART 4: generate ASOS daily netcdf file
# Aim: create a netcdf file of the ASOS station observation data to call and pull from for bias correction
# each station point in netcdf must have their own unique time series and units must match the LOCA2 GCMs

# Step1: load data
# data directory
data_folder = r'Data\ASOS\PreparedData\Daily'
# load ASOS precipitation combined long dataframe
asos_pr_daily_loc = pd.read_csv(path.join(data_folder,'All_ASOS_daily_precipitation_Long.csv'))

# Step 2: set multi-index to use in structuring xarray
asos_pr_daily_loc.set_index(['time', 'lat', 'lon'], inplace=True)

# Step 3: convert to xarray dataset and setup attributes
# convert the dataframe to xarray dataset (ds), those are simply multiple dimension data arrays
asos_daily_pr_ds = xr.Dataset.from_dataframe(asos_pr_daily_loc)
# add attributes to match LOCA2 GCMs and properly describe
asos_daily_pr_ds['pr'].attrs = {'units': 'kg m-2 s-1', 'long_name': 'Precipitation'}
asos_daily_pr_ds['stid'].attrs = {'units': 'N/A', 'long_name': 'Station ID'}
asos_daily_pr_ds.attrs = {'Conventions': 'CF-1.11',
                          'Title': 'ASOS Station pr',
                          'desc': 'ASOS station precipitation (pr) data converted to netcdf format, '
                                  'following CF conventions for units and time'}
# save as netcdf, this should automatically kick it to an xarray dataset
outnc = "ASOS_precipitation_TXYdims.nc"
asos_daily_pr_ds.to_netcdf(path=path.join(data_folder, outnc), format='NETCDF4')

# Finished
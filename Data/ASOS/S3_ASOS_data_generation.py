# -*- coding: utf-8 -*-

"""""
AIM: Prepare precipitation data for each ASOS station
ASOS data downloaded using ASOS_data_collection.py

# 4 Part code for preparing previously extracted ASOS station data for analysis

# PART 1: generate hourly precipitation data
# PART 2: generate event based duration precipitation data
# PART 3: generate yearly precipitation data
# PART 4: generate daily precipitation data
# PART 5: generate monthly precipitation data
# PART 6: generate temporal disaggregation parameters

"""""

import pandas as pd
import json
from os import path
from glob import glob
import csv
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter


# PART 1: Generate hourly precipitation data
# Aim: generate hourly precipitation data for each ASOS station

# set data folder path
asos_folder = r"Data\ASOS\Raw_Precipitation_Data"
# create saving directory if it does not already exist
hourly_dir = r"Data\ASOS\PreparedData\Hourly"
# Check if the directory exists before creating it
if not os.path.exists(hourly_dir):
    os.makedirs(hourly_dir)

# Step 1: list files of raw ASOS precipitation data
asos_files = glob(asos_folder + "/*.csv")
# Check
for name in asos_files:
    print(name)

# Prepare files
for fi in range(len(asos_files)):
    print(fi, "out of", len(asos_files))
    # Step 2: read csvs
    preci = pd.read_csv(asos_files[fi])
    preci.head()

    # Step 3: replace all NA data with 0.0
    preci.iloc[:, 1].fillna(0.0, inplace=True)
    preci.head()

    # Step 4: format time column as datetime format and set as the index
    preci['time'] = pd.to_datetime(preci['time'])
    preci.info()
    preci.set_index('time', inplace=True)

    # Step 5: Because each hour time window recorded the cumulative precipitation,
    # we should use the maximum recorded value within that hour as the hourly precipitation value:
    # Assume that the highest recorded value represents the peak intensity within the hour
    hourly = preci.iloc[:, 0].resample('H').max().fillna(0)

    # Step 6: save files
    site = hourly.columns[0]
    outfile = site + "_hourly_precipitation.csv"
    hourly.to_csv(path.join(hourly_dir, outfile))
    print("finished preparing hourly data for " + site)
    print(outfile)


# PART 2: Generate event based duration precipitation
# Aim: separate the rain events based on duration for each ASOS station

# Step 1: load data from Part 1, hourly precipitation generation
hourly_dir = 'Data/ASOS/PreparedData/Hourly'
##  write folder path to a vector to read at once
hourly_files = glob(hourly_dir + "/*.csv")
# check
len(hourly_files)

# There are many low values as 0.0001, which could be consided as noise.
# Considerting mose values are greater or equal to 0.01. In my following code, I set the rain_threshold to 0.001,
# so that any values smaller than 0.001 will be neglected. -- HN

# Define rain event durations
rain_durations = [1, 2, 3, 4, 6, 12, 24, 48, 72, 96] # in hours

# Define threshold for rainfall (assuming any non-zero precipitation indicates rainfalls)
rain_threshold = 0.0001 # check the distribution for low numbers, set different threshold based on request -- HN

# read csv into a list
hourly_dfs = list()
for fi in range(len(hourly_files)):
    print(fi, "out of", len(hourly_files))
    # read file
    f = hourly_files[fi]
    hourly_precipitation = pd.read_csv(f)

    # Step 2: extract and save the site name
    site = hourly_precipitation.columns[1]

    # Step 3: format time and set as the index
    hourly_precipitation['time'] = pd.to_datetime(hourly_precipitation['time'])
    hourly_precipitation.info()
    hourly_precipitation.set_index('time', inplace=True)
    hourly_precipitation = hourly_precipitation[site]

    # Step 4: Plot
    plt.figure(figsize=(10, 6))
    hourly_precipitation.plot(kind='line', color='blue', linestyle='-', linewidth=2)
    plt.xlabel('Date')
    title = site + 'Hourly Precipitation'
    plt.ylabel('Hourly Precipitation')
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    # check whether directory already exists
    path = "fig/"
    if not os.path.exists(path):
        os.mkdir(path)
        print("Folder %s created!" % path)
    else:
        print("Folder %s already exists" % path)
    plt.savefig(path + title + ".png")

    # Step 5: Separate the rain events based on durations
    # Find the rain events
    ## empty list to hold all stations rain events lists
    rain_events = []
    # Initialized the variables
    rain_start = None
    rain_end = None
    dry_period = 0
    # iterate through the data
    for timestamp, precipitation in hourly_precipitation.items():
        if precipitation > rain_threshold:
            if rain_start is None:
                rain_start = timestamp
            if dry_period >= 6:
                # if there's a dry period of at least 6 hours, consider it a new rain event
                if rain_end:
                    rain_events.append((rain_start, rain_end))
                rain_start = timestamp
                rain_end = None
            dry_period = 0
        else:
            if rain_start is not None:
                dry_period += 1
                if dry_period == 6:
                    # If there's a dry period of at least 6 hours, mark the end of the rain event
                    rain_end = timestamp
                    rain_events.append((rain_start, rain_end))
                    rain_start = None
                    rain_end = None
    # check the number of events
    len(rain_events)
    # adjust the time for rain events
    rain_events_adjusted = [(start, end - pd.Timedelta(hours=5)) for start, end in rain_events]
    # Here I am using a dictionary to save events based on durations,
    # {'1': (start, end, duration, total precipitation)... (start, end, duration,total precipitation), etc;
    #  '2': ...;
    # ....
    # '96': ...}
    filtered_events_by_duration = {}

    for i in range(len(rain_durations)):
        if i == len(rain_durations) - 1:
            upper_bound = rain_durations[i] + 1000
        else:
            upper_bound = rain_durations[i + 1]

        filtered_events = [
            (start, end, int((end - start).total_seconds() / 3600), hourly_precipitation[start:end].sum()) for
            start, end in rain_events_adjusted if
            rain_durations[i] <= (end - start).total_seconds() / 3600 < upper_bound]
        filtered_events_by_duration[str(rain_durations[i])] = filtered_events
        print(f"Rain events lasting {rain_durations[i]}-hr: {len(filtered_events)}")

    # make a new folder path
    path = './PreparedData/Hourly/events_by_durations_' + site
    # check whether directory already exists
    if not os.path.exists(path):
        os.mkdir(path)
        print("Folder %s created!" % path)
    else:
        print("Folder %s already exists" % path)
    # Assuming filtered_events_by_duration is your dictionary
    for key, data in filtered_events_by_duration.items():
        filename = path + f"/{site}_{key}.csv"  # File name based on the key
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Start', 'End', 'Duration', 'Precipitation', 'SiteID'])  # Header row
            for row in data:
                writer.writerow(row)
    print("Events by duration has been completed, and saved to here %s" % path)


# PART 3: Generate yearly precipitation data
# Aim: generate yearly precipitation data for each ASOS station

# Step 1: get data from Part 1, hourly precipitation generation
hourly_dir = 'Data/ASOS/PreparedData/Hourly'
##  write folder path to a vector to read at once
hourly_files = glob(hourly_dir + "/*.csv")
# check
len(hourly_files)
# saving directory for yearly ASOS data
yearly_dir = 'Data/ASOS/PreparedData/Yearly'
# Check if the directory exists before creating it
if not os.path.exists(yearly_dir):
    os.makedirs(yearly_dir)

# prepare csvs in list
for fi in range(len(hourly_files)):
    print(fi, "out of", len(hourly_files))
    # Step 2: read csvs
    preci = pd.read_csv(hourly_files[fi])

    # Step 3: format time and set as the index
    preci['time'] = pd.to_datetime(preci['time'])
    preci.info()
    preci.set_index('time', inplace=True)

    # Step 4: resample the hourly precipitation to get the yearly sum precipitation
    yearly = preci.resample('YE').sum()
    yearly.head()

    # Step 5: save files
    site = yearly.columns[0]
    outfile = site + "_yearly_precipitation.csv"
    yearly.to_csv(path.join(r"PreparedData\Yearly", outfile))

    # Step 6: plot data
    plt.figure(figsize=(10, 6))
    preci.plot(kind='line', color='blue', linestyle='-', linewidth=2)
    plt.xlabel('Year')
    site = preci.columns[0]
    title = site + ' Yearly Precipitation'
    plt.ylabel('Yearly Precipitation (in)')
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()

    # save plot
    # check whether directory already exists
    path = "fig/"
    if not os.path.exists(path):
        os.mkdir(path)
        print("Folder %s created!" % path)
    else:
        print("Folder %s already exists" % path)
    plt.savefig(path + title + ".png")
    print("Finished preparing yearly data and saved to %s" % outfile)


# PART 4: Generate daily precipitation data
# Aim: generate daily precipitation for each ASOS station

# Step 1: load data from Part 1, hourly precipitation generation
hourly_dir = 'Data/ASOS/PreparedData/Hourly'
##  write folder path to a vector to read at once
hourly_files = glob(hourly_dir + "/*.csv")
# check
len(hourly_files)

# saving directory for daily ASOS data
daily_dir = 'Data/ASOS/PreparedData/Daily'
# Check if the directory exists before creating it
if not os.path.exists(daily_dir):
    os.makedirs(daily_dir)

# Prepare daily precipitation dataframes looping through files
for fi in range(len(hourly_files)):
    print(fi, "out of", len(hourly_files))

    # Step 2: read file
    f = hourly_files[fi]
    preci = pd.read_csv(f)

    # Step 3: format time and set as the index
    preci['time'] = pd.to_datetime(preci['time'])
    preci.info()
    preci.set_index('time', inplace=True)

    # Step 4: resample the hourly precipitation to get the yearly sum precipitation
    daily = preci.resample('D').sum()
    print(daily.head())

    # Step 5: reshape
    site = daily.columns[0]
    daily['stid'] = site
    daily = daily.rename(columns = {site:'Precipitation'})

    # Step 6: save the file
    outfile = site + "_daily_precipitation.csv"
    daily.to_csv(path.join(daily_dir, outfile))

    # Step 7: plot data
    plt.figure(figsize=(10, 6))
    preci.plot(kind='line', color='blue', linestyle='-', linewidth=2)
    plt.xlabel('Year')
    site = preci.columns[0]
    title = site + ' Yearly Precipitation'
    plt.ylabel('Yearly Precipitation (in)')
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()

    # save plot
    # check whether directory already exists
    path = "fig/"
    if not os.path.exists(path):
        os.mkdir(path)
        print("Folder %s created!" % path)
    else:
        print("Folder %s already exists" % path)
    plt.savefig(path + title + ".png")
    print(f'{site} daily precipitation is complete and saved')

# PART 5: Generate monthly precipitation data
# Aim: generate monthly precipitation data for each ASOS station
# Step 1: load data from Part 1, hourly precipitation generation
hourly_dir = 'Data/ASOS/PreparedData/Hourly'
##  write folder path to a vector to read at once
hourly_files = glob(data_folder + "/*.csv")
# check
len(hourly_files)
print(hourly_files[:5])

for f in range(len(hourly_files)):
    # Step 1: read csv into a list
    preci = pd.read_csv(hourly_files[f])
    # Step 2: format time and set as the index
    preci['time'] = pd.to_datetime(preci['time'])
    preci = preci.set_index('time')
    # Step 3: resample the hourly precipitation to get the yearly sum precipitation
    monthly = preci.resample('ME').sum()
    # Step 6: save files
    site = monthly.columns[0]
    outfile = site + "_monthly_precipitation.csv"
    monthly.to_csv(path.join(r"PreparedData\Monthly", outfile))
    # Step 5: plot
    monthly.plot(kind='line', color='blue', linestyle='-', linewidth=2, figsize=(10,6))
    plt.xlabel('Date')
    title = site + ' Monthly Precipitation Totals'
    plt.ylabel('Precipitation (in)')
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("fig/" + title + ".png")

# PART 6: Generate ASOS Parameters for Temporal Disaggregation
# AIM: Script to prepare Observation precipitation data for temporal disaggregation using machine learning methods

# Preparing the following parameters:
#     Station Lat, Lon, Elevation
#     Precipitation 1day (P1d)
#     Precipitation -1day (P-1d)
#     Precipitation +1day (P+1d)
#     Precipitation annual total (Pa)
#     Precipitation monthly total (Pm)
#     Precipitation -1month total (P-1m)
#     Precipitation +1month total (P+1m)
# This requires compiling station data from the previously parts of this script

# Step 1: load asos station data that holds lat, lon, and elevation data
# load the saved file to avoid it getting written over
asos_locs_df = pd.read_csv(r'Data\ASOS\TX_ASOS_Stations_CleanClipped.csv')
# get station id's found in asos_final
sites = asos_locs_df['stid'].values.tolist()

# loop through the list sites
for i, stid in enumerate(sites):
    # initial start time
    ist = dt.datetime.now()
    print("Operation start time: ", ist)
    print(i, ":", len(sites)-1)
    print("station ID:", stid)
    # Step 2: Load the station event dataframes for the use of 6-hr and 12-hr disaggregation
    # filepath
    pr_6HRevents_file = ('Data/ASOS/PreparedData/Hourly/' +
                         f'events_by_durations_{stid}/{stid}_6.csv')
    pr_12HRevents_file = ('Data/ASOS/PreparedData/Hourly/' +
                          f'events_by_durations_{stid}/{stid}_12.csv')
    # load files
    pr_6HRevents_df = pd.read_csv(prj_dir + pr_6HRevents_file)
    pr_12HRevents_df = pd.read_csv(prj_dir + pr_12HRevents_file)

    # Step 3: Prepare events dataframe
    # join the two dataframes
    pr_6_12HRevents_df = pd.concat([pr_6HRevents_df, pr_12HRevents_df], ignore_index=True)
    # Correct the SiteID column that was saved empty
    pr_6_12HRevents_df['SiteID'] = stid
    # convert the start and end timestamp into datetime objects
    pr_6_12HRevents_df['Start'] = pd.to_datetime(pr_6_12HRevents_df['Start'])
    pr_6_12HRevents_df['End'] = pd.to_datetime(pr_6_12HRevents_df['End'])
    # Check to see if the start and end dates are on different days
    # diff_days = pr_6_12HRevents_df['End'].dt.day - pr_6_12HRevents_df['Start'].dt.day
    # print('number of events that had different start and end dates ', len(diff_days[diff_days != 0]))
    # keep only the end date as it would be the day to merge with the others
    pr_6_12HRevents_to_merge_df = pr_6_12HRevents_df.drop(columns='Start').copy()
    # Create dataframe
    pr_6_12HRevents_enddates_df = pd.DataFrame({'End': pr_6_12HRevents_to_merge_df['End'].dt.date,
                                                'Year': pr_6_12HRevents_to_merge_df['End'].dt.year,
                                                'Month': pr_6_12HRevents_to_merge_df['End'].dt.month,
                                                'Day': pr_6_12HRevents_to_merge_df['End'].dt.day})
    # Step 4: Load the observational yearly, monthly, and daily datasets
    # filepaths
    monthly_obs_file = ('Data/ASOS/PreparedData/Monthly/' +
                        f'{stid}_monthly_precipitation.csv')
    yearly_obs_file = ('Data/ASOS/PreparedData/Yearly/' +
                       f'{stid}_yearly_precipitation.csv')
    daily_obs_file = ('Data/ASOS/PreparedData/Daily/' +
                      f'{stid}_daily_precipitation.csv')
    # load files
    monthly_obs_df = pd.read_csv(prj_dir + monthly_obs_file)
    yearly_obs_df = pd.read_csv(prj_dir + yearly_obs_file)
    daily_obs_df = pd.read_csv(prj_dir + daily_obs_file)

    # Step 5: Prepare the observational yearly, monthly, and daily datasets

    # Prepare the dataframes for datetime objects
    monthly_obs_df['time'] = pd.to_datetime(monthly_obs_df['time'])
    yearly_obs_df['time'] = pd.to_datetime(yearly_obs_df['time'])
    daily_obs_df['time'] = pd.to_datetime(daily_obs_df['time'])

    # Prepare the monthly rainfall data
    monthly_obs_df['Year'] = monthly_obs_df['time'].dt.year
    monthly_obs_df['Month'] = monthly_obs_df['time'].dt.month
    # Create the events monthly totals dataframe (Pm)
    # join by year and month
    pr_events_monthlyTotal_df = pd.merge(pr_6_12HRevents_enddates_df, monthly_obs_df, on=['Year', 'Month'], how='left')
    # clean the table to prepare for final merge
    pr_events_monthlyTotal = pr_events_monthlyTotal_df[stid].rename('pr_monthly_T')

    # Prepare the yearly rainfall data
    yearly_obs_df['Year'] = yearly_obs_df['time'].dt.year
    # Create the events yearly totals dataframe (Pa)
    # join by year
    pr_events_yearlyTotal_df = pd.merge(pr_6_12HRevents_enddates_df, yearly_obs_df, on=['Year'], how='left')
    # clean the table to prepare for final merge
    pr_events_yearlyTotal = pr_events_yearlyTotal_df[stid].rename('pr_annual_T')

    # Prepare the daily rainfall data
    # Create the events daily totals dataframe (P1d)
    pr_6_12HRevents_enddates_df['End'] = pd.to_datetime(pr_6_12HRevents_enddates_df['End'])
    # join by date
    pr_events_1dayTotal_df = pd.merge(pr_6_12HRevents_enddates_df, daily_obs_df,
                                      how="left", left_on="End", right_on="time")
    # clean the table to prepare for final merge
    pr_events_1dayTotal_df = pr_events_1dayTotal_df.rename(columns={'Precipitation': 'pr_1day_T'})
    # pr_events_1dayTotal_df = pr_events_1dayTotal_df.loc[:, ['time', 'pr_1day_T']]
    pr_events_1dayTotal = pr_events_1dayTotal_df['pr_1day_T']
    # Create the events +1 day totals dataframe (P+1d)
    # create a time offset
    offset = pd.Timedelta(1, 'D')
    # get the +1 day dates
    pr_events_add1dates = pr_6_12HRevents_enddates_df['End'] + offset
    # join by date
    pr_events_add1dayTotal_df = pd.merge(pr_events_add1dates, daily_obs_df,
                                         how="left", left_on="End", right_on="time")
    # pull the value and rename to prepare for final merge
    pr_events_add1dayTotal_df = pr_events_add1dayTotal_df.rename(columns={'Precipitation': 'pr_+1day_T'})
    # pr_events_add1dayTotal_df = pr_events_add1dayTotal_df.loc[:, ['time', 'pr_+1day_T']]
    pr_events_add1dayTotal = pr_events_add1dayTotal_df['pr_+1day_T']
    # Create the events +1 day totals dataframe (P-1d)
    # get the -1 day dates
    pr_events_minus1dates = pr_6_12HRevents_enddates_df['End'] - offset
    # join by date
    pr_events_minus1dayTotal_df = pd.merge(pr_events_minus1dates, daily_obs_df,
                                           how="left", left_on="End", right_on="time")
    # clean the table to prepare for final merge
    pr_events_minus1dayTotal_df = pr_events_minus1dayTotal_df.rename(columns={'Precipitation': 'pr_-1day_T'})
    # pr_events_minus1dayTotal_df = pr_events_minus1dayTotal_df.loc[:, ['time', 'pr_-1day_T']]
    pr_events_minus1dayTotal = pr_events_minus1dayTotal_df['pr_-1day_T']
    # Create the events +1 Month dataframe (P+1m)
    # get the +1 month dates
    events_add1month = pr_6_12HRevents_enddates_df['End'] + pd.offsets.DateOffset(months=1)
    # create a new dataframe to hold the information
    add1month_to_merge_df = pd.DataFrame({'Year': events_add1month.dt.year,
                                          'Month': events_add1month.dt.month})
    # merge with the Monthly Totals dataframe to get the values for column
    pr_events_add1month_df = pd.merge(add1month_to_merge_df, monthly_obs_df,
                                      how='left', on=['Year', 'Month'])
    # drop the extra columns and rename
    pr_events_add1month = pr_events_add1month_df[stid].rename('pr_+1monthly_T')
    # Create the events -1 Month dataframe (P-1m)
    # get the -1 month dates
    events_minus1month = pr_6_12HRevents_enddates_df['End'] - pd.offsets.DateOffset(months=1)
    # create a new dataframe to hold the information
    minus1month_to_merge_df = pd.DataFrame({'Year': events_minus1month.dt.year,
                                            'Month': events_minus1month.dt.month})
    # merge with the Monthly Totals dataframe to get the values for column
    pr_events_minus1month_df = pd.merge(minus1month_to_merge_df, monthly_obs_df,
                                        how='left', on=['Year', 'Month'])
    # drop the extra columns and rename
    pr_events_minus1month = pr_events_minus1month_df[stid].rename('pr_-1monthly_T')

    # Step 4: Prepare the training, validation, testing dataframe
    # rename events precip column adding '_y' to easily identify the predictand
    pr_6_12HRevents_to_merge_df = pr_6_12HRevents_to_merge_df.rename(columns={'Precipitation': 'Precipitation_y'})
    # split the time and date of the end column
    pr_6_12HRevents_to_merge_df['End_Date'] = pr_6_12HRevents_to_merge_df['End'].dt.date
    pr_6_12HRevents_to_merge_df['End_Time'] = pr_6_12HRevents_to_merge_df['End'].dt.time
    # create a day of year column
    pr_6_12HRevents_to_merge_df['DOY'] = pr_6_12HRevents_to_merge_df['End'].dt.dayofyear
    # get the station details
    stid_locs = asos_locs_df.loc[asos_locs_df['stid'] == stid, ['stid', 'lat', 'lon', 'elev']]
    # Add station details to final dataframe by station ID
    pr_6_12HRevents_locs_df = pd.merge(pr_6_12HRevents_to_merge_df, stid_locs, left_on="SiteID",
                                       right_on="stid", how="left")
    # Add monthly totals to final dataframe by index
    pr_6_12HRevents_locs_m_df = pd.concat([pr_6_12HRevents_locs_df, pr_events_monthlyTotal], axis=1)
    # Add annual totals to final dataframe by index
    pr_6_12HRevents_locs_ma_df = pd.concat([pr_6_12HRevents_locs_m_df, pr_events_yearlyTotal], axis=1)
    # Add 1day totals to final dataframe by index
    pr_6_12HRevents_locs_mad_df = pd.concat([pr_6_12HRevents_locs_ma_df, pr_events_1dayTotal], axis=1)
    # Add +1day totals to final dataframe by index
    pr_6_12HRevents_locs_madp_df = pd.concat([pr_6_12HRevents_locs_mad_df, pr_events_add1dayTotal], axis=1)
    # Add -1day totals to final dataframe by index
    pr_6_12HRevents_locs_madpm_df = pd.concat([pr_6_12HRevents_locs_madp_df, pr_events_minus1dayTotal], axis=1)
    # Add +1month totals to final dataframe by index
    pr_6_12HRevents_locs_madp2m_df = pd.concat([pr_6_12HRevents_locs_madpm_df, pr_events_add1month], axis=1)
    # Add -1month totals to final dataframe by index
    pr_6_12HRevents_locs_madp2m2_df = pd.concat([pr_6_12HRevents_locs_madp2m_df, pr_events_minus1month], axis=1)
    # check the columns exist
    print(pr_6_12HRevents_locs_madp2m2_df.columns)
    # drop excess station ID column
    pr_6_12HRevents_locs_final_df = pr_6_12HRevents_locs_madp2m2_df.drop(columns=["SiteID", "End"])
    # round to only two digit integers for most columns
    cols_round = ['elev', 'Precipitation_y', 'pr_monthly_T', 'pr_annual_T', 'pr_1day_T', 'pr_+1day_T', 'pr_-1day_T',
                  'pr_+1monthly_T', 'pr_-1monthly_T']
    pr_6_12HRevents_locs_final_df[cols_round] = pr_6_12HRevents_locs_final_df[cols_round].round(2)

    # arrange columns by certain order
    col_order = ['stid', 'End_Date', 'End_Time', 'Duration', 'Precipitation_y', 'DOY', 'lat', 'lon', 'elev',
                 'pr_monthly_T', 'pr_annual_T', 'pr_1day_T', 'pr_+1day_T', 'pr_-1day_T',
                 'pr_+1monthly_T', 'pr_-1monthly_T']
    pr_6_12HRevents_locs_final_df = pr_6_12HRevents_locs_final_df[col_order]
    # create the 6-HR events dataframe
    pr_6 = pr_6_12HRevents_locs_final_df[pr_6_12HRevents_locs_final_df['Duration'] == 6]
    # create the 12-HR events dataframe
    pr_12 = pr_6_12HRevents_locs_final_df[pr_6_12HRevents_locs_final_df['Duration'] == 12]
    # save the dataframes to where it will be used
    pr_6_12HRevents_locs_final_df.to_csv(f'Temporal disaggregation/{stid}_6_12_ASOSparameters.csv',
                                         index=False)
    pr_6.to_csv(f'Temporal disaggregation/{stid}_6_ASOSparameters.csv', index=False)
    pr_12.to_csv(f'Temporal disaggregation/{stid}_12_ASOSparameters.csv', index=False)
    del monthly_obs_df, yearly_obs_df, daily_obs_df, pr_events_monthlyTotal_df, pr_events_monthlyTotal, (
        pr_events_yearlyTotal_df), pr_events_yearlyTotal, pr_events_1dayTotal_df, pr_events_1dayTotal, (
        pr_events_add1dates), pr_events_add1dayTotal_df, pr_events_add1dayTotal, pr_events_minus1dates, (
        pr_events_minus1dayTotal_df), pr_events_minus1dayTotal, events_add1month, add1month_to_merge_df, (
        pr_events_add1month_df), pr_events_add1month, events_minus1month, minus1month_to_merge_df, (
        pr_events_minus1month_df), pr_events_minus1month
    del pr_12, pr_6, pr_6_12HRevents_locs_final_df, pr_6_12HRevents_locs_m_df, pr_6_12HRevents_locs_df, (
        pr_6_12HRevents_locs_ma_df), pr_6_12HRevents_locs_mad_df, pr_6_12HRevents_locs_madp_df, (
        pr_6_12HRevents_locs_madpm_df), pr_6_12HRevents_locs_madp2m_df, pr_6_12HRevents_locs_madp2m2_df, (
        stid_locs), pr_6_12HRevents_to_merge_df, pr_6_12HRevents_df
    gc.collect()
    # initial start time
    endt = dt.datetime.now()
    print("length of operation: ", endt-ist)

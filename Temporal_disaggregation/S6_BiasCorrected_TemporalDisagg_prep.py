# -*- coding: utf-8 -*-
"""""
AIM: Script to prepare Bias Corrected climate data for temporal disaggregation using machine learning

Preparing the following parameters:
    Station Lat, Lon, Elevation
    Precipitation 1day (P1d)
    Precipitation -1day (P-1d)
    Precipitation +1day (P+1d)
    Precipitation annual total (Pa)
    Precipitation monthly total (Pm)
    Precipitation -1month total (P-1m)
    Precipitation +1month total (P+1m)


1 Part code for preparing bias corrected for machine learning models prepare bias corrected GCMs data by generating 
new parameters.

"""""
import gc
import pandas as pd
import glob as glob
from os import getcwd
import datetime as dt
from matplotlib import pyplot as plt

# Function to clean filepath list
# Removes the excess csv files that are either not a completed model set, the final list of sites, bias correction
# methods that are not completed
# requires two inputs:
# list_filepaths: list of file pathways to be cleaned
# siteIDs: list of station IDs
def clean_bc_filepath_list(list_filepaths, siteIDs):
    # Removes delta method
    cleaned_list = [f for f in list_filepaths if not f.endswith('delta_method.csv')]
    # Removes alternative last step applied to methods that used timesteps
    cleaned_list = [f for f in cleaned_list if not '_altLastTStep' in f]
    # Remove detrended quantile mapping method
    cleaned_list = [f for f in cleaned_list if not 'detrended' in f]
    # Remove the GCM INM-CM5-0
    cleaned_list = [f for f in cleaned_list if not 'INM-CM5-0' in f]
    # Remove stations that weren't in the final list
    cleaned_list = [i for e in siteIDs for i in cleaned_list if e in i]
    return cleaned_list

# Function to create the file name
def create_bc_ml_filename(filepath):
    # get string index to use for unique identifiers
    str_i = filepath.find(
        '5\\') + 2  # '5\\' is unique string found at the end of all files before naming convention
    # name string to use add the application reference and file type
    return filepath[str_i:-4] + '_MLinput.csv'

# Function to get the station ID from the file name
def get_stid(filepath):
    # get string index to use for unique identifiers
    str_i = filepath.find(
        '5\\') + 2  # '5\\' is unique string found at the end of all files before naming convention
    # name string to use add the application reference and file type
    return filepath[str_i:str_i+3]

# initial start time
ist = dt.datetime.now()
print("Initial script start time: ", ist)

# Step 1: Load data and setup workspace

# get current working directroy
cwd_path = getcwd()

# load ASOS station lat, lon, elev description dataframe
asos_locs_df = pd.read_csv(cwd_path + '/Data/ASOS/TX_ASOS_Final.csv')

# list station id's found in asos_final
sites = asos_locs_df['stid'].values.tolist()

# list bias corrected filepaths
ssp245_files = glob.glob(cwd_path + "BiasCorrection/BiasCorrected_Data/*_biascorr/SSP245/*.csv", recursive=True)
ssp585_files = glob.glob(cwd_path + "BiasCorrection/BiasCorrected_Data/*_biascorr/SSP585/*.csv", recursive=True)
all_files = ssp245_files + ssp585_files
cleaned_files = clean_bc_filepath_list(list_filepaths=all_files, siteIDs=sites)
len(cleaned_files) # this should be divisible by the number of station IDs

# create a list of stations that may have potential errors
double_check_stids = []
# Loop through stations to prepare data
for i, file in enumerate(cleaned_files):

    # initial start time
    ist = dt.datetime.now()
    print("Operation start time: ", ist)
    print(i, ":", len(cleaned_files)-1)

    # Step 2: create saving filepaths
    # get the station ID for this file
    stid = get_stid(file)
    print(stid)
    # create saving filename
    fname = create_bc_ml_filename(file)
    fpath = 'Temporal disaggregation/' + fname  # adds directory pathway to save to
    # bypass if the file already exists
    if os.path.exists(fpath):
        print(fname, "already exists; SKIP")
    else:
        # Step 3: Load bias corrected precip data files and prepare for generating parameters by isolating to events

        # load file
        bias_corr_df = pd.read_csv(file)
        # create datetime objects
        bias_corr_df['time'] = pd.to_datetime(bias_corr_df['time'])
        # create datetime index
        bias_corr_df = bias_corr_df.set_index('time')
        # first check to make sure there are no negative values
        min_val = bias_corr_df.min().values
        # if there are negative values change them to 0
        if min_val < 0:
            print("negative minimum precipitation ", min_val)
            # Replace negative values with 0
            bias_corr_df = bias_corr_df.map(lambda x: max(x, 0))
        else:
            print("non-negative minimum precipitation")

        # Prepare the daily rainfall data as events
        # set minimum threshold based on event durations data generation
        rain_threshold = 0.01
        # remove all 0 rainfall days
        bias_corr_daily_df = bias_corr_df[bias_corr_df['pr'] > rain_threshold]

        # create an alert to identify if the first/end date are the first/last date of whole timeseries
        if bias_corr_daily_df.index.min() == pd.Timestamp(dt.datetime(2015, 1, 1)):
            print("potential error first day of rain above threshold is first day\nStation", stid)
            double_check_stids.append({'first', stid})
            # drop this row
            bias_corr_daily_df = bias_corr_daily_df.drop(index=bias_corr_daily_df.index.min())
        if bias_corr_daily_df.index.max() == pd.Timestamp(dt.datetime(2100, 12, 31)):
            print("potential error last day of rain above threshold is last day\nStation", stid)
            double_check_stids.append({'last', stid})
            # drop this row
            bias_corr_daily_df = bias_corr_daily_df.drop(index=bias_corr_daily_df.index.max())

        # Create the events daily totals dataframe (P1d), base to build dataframe from
        # clean the table to prepare for final merge
        bias_corr_daily_df = bias_corr_daily_df.rename(columns={'pr': 'pr_1day_T'})

        # Step 4: Prepare climate parameters

        # Create the events monthly totals dataframe (Pm)
        # calculate the monthly totals
        pr_monthlyTotal_df = bias_corr_df.resample('ME').sum()
        # prepare for future join by adding year column
        pr_monthlyTotal_df['Month'] = pr_monthlyTotal_df.index.month
        # prepare for future join by adding year column
        pr_monthlyTotal_df['Year'] = pr_monthlyTotal_df.index.year
        # rename
        pr_monthlyTotal_df = pr_monthlyTotal_df.rename(columns={'pr': 'pr_monthly_T'})
        # reset and drop index
        pr_monthlyTotal_df = pr_monthlyTotal_df.reset_index()
        pr_monthlyTotal_df = pr_monthlyTotal_df.drop(columns='time')

        # Create the events yearly totals dataframe (Pa)
        # Prepare the yearly rainfall data
        pr_yearlyTotal_df = bias_corr_df.resample('YE').sum()
        # prepare for future join by adding year column
        pr_yearlyTotal_df['Year'] = pr_yearlyTotal_df.index.year
        # rename
        pr_yearlyTotal_df = pr_yearlyTotal_df.rename(columns={'pr': 'pr_annual_T'})
        # reset and drop date
        pr_yearlyTotal_df = pr_yearlyTotal_df.reset_index()
        pr_yearlyTotal_df = pr_yearlyTotal_df.drop(columns='time')

        # Create the events +1 day totals dataframe (P+1d)
        # create a time offset
        offset = pd.Timedelta(1, 'D')
        # get the +1 day dates
        pr_events_add1dates = bias_corr_daily_df.index + offset
        # keep by date
        pr_events_add1dayTotal_df = bias_corr_df[bias_corr_df.index.isin(pr_events_add1dates)]
        # pull the value and rename to prepare for final merge
        pr_events_add1dayTotal = pr_events_add1dayTotal_df.rename(columns={'pr': 'pr_+1day_T'})
        # reset and remove date index
        pr_events_add1dayTotal = pr_events_add1dayTotal.reset_index()
        pr_events_add1dayTotal = pr_events_add1dayTotal.drop(columns='time')
        # Create the events +1 day totals dataframe (P-1d)
        # get the -1 day dates
        pr_events_minus1dates = bias_corr_daily_df.index - offset
        # join by date
        pr_events_minus1dayTotal_df = bias_corr_df[bias_corr_df.index.isin(pr_events_minus1dates)]
        # clean the table to prepare for final merge
        pr_events_minus1dayTotal = pr_events_minus1dayTotal_df.rename(columns={'pr': 'pr_-1day_T'})
        # reset and remove date index
        pr_events_minus1dayTotal = pr_events_minus1dayTotal.reset_index()
        pr_events_minus1dayTotal = pr_events_minus1dayTotal.drop(columns='time')

        # Create the events +1 Month dataframe (P+1m)
        # get the +1 month dates
        events_add1month = bias_corr_daily_df.index + pd.offsets.DateOffset(months=1)
        # create a new dataframe to hold the information
        add1month_to_merge_df = pd.DataFrame({'Year': events_add1month.year,
                                              'Month': events_add1month.month})
        # merge with the Monthly Totals dataframe to get the values for column
        pr_events_add1month_df = pd.merge(add1month_to_merge_df, pr_monthlyTotal_df, how='left', on=['Year', 'Month'])
        # drop the extra columns and rename
        pr_events_add1month = pr_events_add1month_df['pr_monthly_T'].rename('pr_+1monthly_T')

        # Create the events -1 Month dataframe (P-1m)
        # get the -1 month dates
        events_minus1month = bias_corr_daily_df.index - pd.offsets.DateOffset(months=1)
        # create a new dataframe to hold the information
        minus1month_to_merge_df = pd.DataFrame({'Year': events_minus1month.year,
                                                'Month': events_minus1month.month})
        # merge with the Monthly Totals dataframe to get the values for column
        pr_events_minus1month_df = pd.merge(minus1month_to_merge_df, pr_monthlyTotal_df, how='left', on=['Year', 'Month'])
        # drop the extra columns and rename
        pr_events_minus1month = pr_events_minus1month_df['pr_monthly_T'].rename('pr_-1monthly_T')

        # Step 5: Prepare the data to create the bias correction temporal disaggregation dataset

        # Add year and month columns
        bias_corr_daily_df['Date'] = bias_corr_daily_df.index.date
        bias_corr_daily_df['Month'] = bias_corr_daily_df.index.month
        bias_corr_daily_df['Year'] = bias_corr_daily_df.index.year
        bias_corr_daily_df['DOY'] = bias_corr_daily_df.index.dayofyear
        # reset the index
        bias_corr_daily_df = bias_corr_daily_df.reset_index()
        # Add station ID used to merge on
        bias_corr_daily_df['stid'] = stid
        # get the station details
        stid_locs = asos_locs_df.loc[asos_locs_df['stid'] == stid, ['stid', 'lat', 'lon', 'elev']]

        # Step 6: Join parameters to base dataframe

        # Add station details to final dataframe by station ID
        bias_corr_daily_locs_df = pd.merge(bias_corr_daily_df, stid_locs,
                                           on="stid", how="left")
        # Add monthly totals to final dataframe by year and month
        bias_corr_events_locs_m_df = pd.merge(bias_corr_daily_locs_df, pr_monthlyTotal_df, on=['Year', 'Month'])
        # Add annual totals to final dataframe by year
        bias_corr_events_locs_ma_df = pd.merge(bias_corr_events_locs_m_df, pr_yearlyTotal_df, on='Year')
        # Add +1day totals to final dataframe by index
        bias_corr_events_locs_map_df = pd.concat([bias_corr_events_locs_ma_df, pr_events_add1dayTotal], axis=1)
        # Add -1day totals to final dataframe by index
        bias_corr_events_locs_madpm_df = pd.concat([bias_corr_events_locs_map_df, pr_events_minus1dayTotal], axis=1)
        # Add +1month totals to final dataframe by index
        bias_corr_events_locs_madp2m_df = pd.concat([bias_corr_events_locs_madpm_df, pr_events_add1month], axis=1)
        # Add -1month totals to final dataframe by index
        bias_corr_events_locs_madp2m2_df = pd.concat([bias_corr_events_locs_madp2m_df, pr_events_minus1month], axis=1)
        # drop the columns no longer needed
        bias_corr_events_locs_final_df = bias_corr_events_locs_madp2m2_df.drop(columns=['Year', 'Month', 'time'])
        # round to only two digit integers for most columns
        cols_round = ['elev', 'pr_monthly_T', 'pr_annual_T', 'pr_1day_T', 'pr_+1day_T', 'pr_-1day_T', 'pr_+1monthly_T',
                      'pr_-1monthly_T']
        bias_corr_events_locs_final_df[cols_round] = bias_corr_events_locs_final_df[cols_round].round(2)
        # organize column order
        col_order = ['stid', 'Date', 'DOY', 'lat', 'lon', 'elev',
                     'pr_monthly_T', 'pr_annual_T', 'pr_1day_T', 'pr_+1day_T', 'pr_-1day_T',
                     'pr_+1monthly_T', 'pr_-1monthly_T']
        bias_corr_events_locs_final_df = bias_corr_events_locs_final_df[col_order]
        # save the dataframe
        bias_corr_events_locs_final_df.to_csv(fpath, index=False)

        # cleanup memory
        del pr_monthlyTotal_df, pr_yearlyTotal_df, pr_events_add1dates, pr_events_add1dayTotal_df, (
            pr_events_add1dayTotal), (
            pr_events_minus1dates), pr_events_minus1dayTotal_df, pr_events_minus1dayTotal, events_add1month, (
            add1month_to_merge_df), pr_events_add1month_df, pr_events_add1month, events_minus1month, (
            minus1month_to_merge_df), pr_events_minus1month_df, pr_events_minus1month
        del bias_corr_events_locs_final_df, bias_corr_df, bias_corr_daily_df, bias_corr_events_locs_map_df, (
            bias_corr_daily_locs_df), bias_corr_events_locs_m_df, bias_corr_events_locs_ma_df, (
            bias_corr_events_locs_madpm_df), bias_corr_events_locs_madp2m_df, bias_corr_events_locs_madp2m2_df, (
            stid_locs)
        gc.collect()
    # initial start time
    endt = dt.datetime.now()
    print("length of operation: ", endt-ist)

# check the list of stations that may have errors related to start/end dates have precipitation
print(double_check_stids)

# Finished
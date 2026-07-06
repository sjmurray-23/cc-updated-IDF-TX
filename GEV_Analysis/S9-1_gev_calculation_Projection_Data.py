# -*- coding: utf-8 -*-
"""
This script is writen to compute the GEV and gumble ditribution for  extreme precipitation projections

3 Part script that performs extreme value distribution analysis on 24-HR, 12-HR, and 6-HR
    PART 1: 24-HR bias corrected GCM projected data
    PART 2: 12-HR bias corrected GCM projected data
    PART 3: 6-HR bias corrected GCM projected data
"""
import os
import pandas as pd
import numpy as np
import xarray as xr
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import genextreme as gev
from scipy.stats import gumbel_r

from glob import glob


###Defining function to compute empirical return level
def empirical_return_level(data):
    df = pd.DataFrame(index=np.arange(data.size))
    # sort the data
    df["sorted"] = np.sort(data)[::-1]
    # rank via scipy instead to deal with duplicate values
    df["ranks_sp"] = np.sort(stats.rankdata(-data))
    # find exceedence probability
    n = data.size
    df["exceedance"] = df["ranks_sp"] / (n + 1)
    # find return period
    df["period"] = 1 / df["exceedance"]

    df = df[::-1]

    out = xr.DataArray(
        dims=["period"],
        coords={"period": df["period"]},
        data=df["sorted"],
        name="level",
    )
    return out


#### Direction file. Populate with the directory address
# get current working directory
cwd_path = os.getcwd()
dir = os.path.join(cwd_path, "GEV_Analysis")
dir_24master_biasCorrected = os.path.join(cwd_path,
                                          "GEV_Analysis/Station_24HR_Maximums/Station_biascorr_24HRmax")
dir_12master_biasCorrected = os.path.join(cwd_path,
                                          "GEV_Analysis/Station_12HR_Maximums/Station_biascorr_12HRmax")
dir_6master_biasCorrected = os.path.join(cwd_path, "GEV_Analysis/Station_6HR_Maximums/Station_biascorr_6HRmax")

##listing scenarios anb models
Models = ['CanESM5', 'GFDL-ESM4', 'BCC-CSM2-MR', 'EC-Earth3', 'CNRM-CM6-1', 'ACCESS-CM2', 'CNRM-ESM2-1']
Methods = ["addition linear scaling", "addition quantile delta mapping", "addition quantile mapping", \
           "multiplication linear scaling", "multiplication quantile delta mapping", "multiplication quantile mapping"]
Scenarios = ["SSP245", "SSP585"]

## reading starting date and the number of years from the stations
dir_master_ASOS = os.path.join(cwd_path, '/Data/ASOS/')
file_path_ASOS = os.path.join(dir_master_ASOS, 'TX_ASOS_Final.csv')

df_ASOS_stations = pd.read_csv(file_path_ASOS)
merged_StationList = df_ASOS_stations['stid'].tolist()

#### Part 1: 24-HR bias corrected GCM projected data

sumArraySeries = []
for station in merged_StationList:
    file_paths = glob(os.path.join(dir_24master_biasCorrected, station + "_bias_corr_24HRmax_df.csv"))

    for file_path in file_paths:
        df_station = pd.read_csv(file_path)
        # Extract the file name from the path
        file_name = os.path.basename(file_path).split('.')[0]

        # Define the start date and the number of days
        start_date = "2015-01-01"
        num_days = 31411
        num_years = 86
        start_year = 2015

        for Scenario in Scenarios:
            Dict_scenario = {}
            for model in Models:
                GEV_out_df = pd.DataFrame(columns=['0.9', '0.99'], index=Methods)
                OutfileName_GEV = station + "_GEV_" + model + "_" + Scenario

                Gumble_out_df = pd.DataFrame(columns=['0.9', '0.99'], index=Methods)
                OutfileName_Gumble = station + "_Gumbel_" + model + "_" + Scenario
                for method in Methods:
                    df_filtered = df_station.loc[(df_station['model'] == model) & \
                                                 (df_station['methods'] == method) & \
                                                 (df_station['scenario'] == Scenario), \
                        ['year', 'pr', "stid"]]

                    # Generate the date range
                    time_index = pd.date_range(start=start_year, periods=num_years, freq="Y")

                    # Create the DataFrame with this index
                    df = pd.DataFrame(index=time_index)
                    df.index.name = "year"
                    df['pr'] = None
                    df['pr'] = df_filtered['pr'].values

                    data = df
                    precipitation = data.pr

                    shape, loc, scale = gev.fit(precipitation.values, 0)

                    if np.all(precipitation.values == precipitation.values[0]):
                        continue
                    loc_g, scale_g = gumbel_r.fit(precipitation.values)

                    # create vector of years
                    years = np.arange(1.1, 100, 0.1)

                    output_dir = dir
                    Output_filename = station + "_" + model + "_" + method + "_" + Scenario

                    # calcuting the probablities of GEV and Gumble distributions

                    vals = gev.ppf([0.9, 0.99], shape, loc=loc)
                    GEV_out_df.loc[method] = vals

                    vals = gumbel_r.ppf([0.9, 0.99], loc=loc_g)
                    Gumble_out_df.loc[method] = vals

                pieces = {"GEV": GEV_out_df, "Gumble": Gumble_out_df}
                result_scenario_model = pd.concat(pieces)

                Dict_scenario[model] = result_scenario_model
            Dict_scenario_df = pd.concat(Dict_scenario)

            ### From this point on it is added to create the summary of the GEV and gumble distributions
            print("Station and Scenarios are :", station, Scenario)
            print()
            print(" 90th percentile mean, 5th, and 95th,  ")
            print(Dict_scenario_df["0.9"].mean(), np.percentile(Dict_scenario_df["0.9"], 5),
                  np.percentile(Dict_scenario_df["0.9"], 95))
            SumArray = [station, Scenario, "0.9", Dict_scenario_df["0.9"].mean(),
                        np.percentile(Dict_scenario_df["0.9"], 5), np.percentile(Dict_scenario_df["0.9"], 95)]
            sumArraySeries.append(SumArray)
            print()
            print("Station and Scenarios are :", station, Scenario)
            print(" 99th percentile mean, 5th, and 95th,  ")
            print(Dict_scenario_df["0.99"].mean(), np.percentile(Dict_scenario_df["0.99"], 5),
                  np.percentile(Dict_scenario_df["0.99"], 95))
            SumArray = [station, Scenario, "0.99", Dict_scenario_df["0.99"].mean(),
                        np.percentile(Dict_scenario_df["0.99"], 5), np.percentile(Dict_scenario_df["0.99"], 95)]
            sumArraySeries.append(SumArray)
            print()

    Summary_df = pd.DataFrame(data=np.array(sumArraySeries),
                              columns=['Station', 'Scenario', 'return_period', 'mean', '5th_percentile',
                                       '95th_percentile'])

separator = ", "

Summary_df.to_csv('All_stations_GCM_IDF_24hr_summary.csv')

#### PART 2: 12-HR bias corrected GCM projected data

# 12-HR temporally disaggregated bias corrected data
sumArraySeries = []
dict_not_enough_years = {}
all_missing = 0

for si in range(len(merged_StationList)):
    print(si, ":", len(merged_StationList))
    station = merged_StationList[si]
    file_paths = glob(os.path.join(dir_12master_biasCorrected, station + "_bias_corr_12HRmax_df.csv"))

    if len(file_paths) == 0:
        pass
    else:
        # Import the bias corrected rainfall data:
        for file_path in file_paths:

            df_station = pd.read_csv(file_path)
            # Extract the file name from the path
            file_name = os.path.basename(file_path).split('.')[0]

            # Define the start date and the number of days
            start_date = "2015-01-01"
            num_days = 31411
            num_years = 86
            start_year = "2015"

            for Scenario in Scenarios:
                Dict_scenario = {}
                for model in Models:
                    GEV_out_df = pd.DataFrame(columns=['0.9', '0.99'], index=Methods)
                    OutfileName_GEV = station + "_GEV_" + model + "_" + Scenario

                    Gumble_out_df = pd.DataFrame(columns=['0.9', '0.99'], index=Methods)
                    OutfileName_Gumble = station + "_Gumbel_" + model + "_" + Scenario
                    for method in Methods:
                        df_filtered = df_station.loc[(df_station['model'] == model) & \
                                                     (df_station['methods'] == method) & \
                                                     (df_station['scenario'] == Scenario), \
                            ['year', 'Predicted_12hr_pr', "stid"]]

                        # Generate the date range
                        time_index = pd.date_range(start=start_year, freq="YS", periods=num_years)

                        # Create the DataFrame with this index
                        df = pd.DataFrame(index=time_index)
                        df.index.name = "year"
                        df['Predicted_12hr_pr'] = None
                        try:
                            df['Predicted_12hr_pr'] = df_filtered['Predicted_12hr_pr'].values
                        except ValueError:
                            print("Length of values did not match, dataset years < 86")
                        if len(df_filtered['Predicted_12hr_pr'].values) < 86:
                            print("correcting the mismatching length value issue by filling forward")
                            df_filtered['year'] = pd.to_datetime(df_filtered['year'], format='%Y')
                            df_filtered_fullrange = df_filtered.set_index('year').reindex(time_index)
                            # sum the number of missing years
                            n_missing = 86 - len(df_filtered['Predicted_12hr_pr'].values)
                            all_missing = all_missing + n_missing
                            print("Number of missing years measurements", all_missing)
                            not_enough_measures = {
                                all_missing: {'stid': station,
                                              'number_missing_years': n_missing,
                                              'missing_year': df_filtered_fullrange[df_filtered_fullrange[
                                                  'Predicted_12hr_pr'].isnull()].index.year.values,
                                              'scenario': Scenario,
                                              'model': model,
                                              'methods': method}
                            }
                            print(not_enough_measures)
                            dict_not_enough_years.update(not_enough_measures)

                            # fill forward the missing information
                            df_filtered_ff = df_filtered_fullrange.ffill()
                            df['Predicted_12hr_pr'] = df_filtered_ff['Predicted_12hr_pr'].values
                        else:
                            df['Predicted_12hr_pr'] = df_filtered['Predicted_12hr_pr'].values

                        #### this block is for computing the max values
                        data = df
                        precipitation = data.Predicted_12hr_pr

                        shape, loc, scale = gev.fit(precipitation.values, 0)

                        if np.all(precipitation.values == precipitation.values[0]):
                            continue
                        loc_g, scale_g = gumbel_r.fit(precipitation.values)

                        # create vector of years
                        years = np.arange(1.1, 100, 0.1)

                        output_dir = dir
                        Output_filename = station + "_" + model + "_" + method + "_" + Scenario

                        # calcuting the probablities of GEV and Gumble distributions

                        vals = gev.ppf([0.9, 0.99], shape, loc=loc)
                        GEV_out_df.loc[method] = vals

                        vals = gumbel_r.ppf([0.9, 0.99], loc=loc_g)
                        Gumble_out_df.loc[method] = vals

                    # these two lines are comments for memory save

                    # GEV_out_df.to_csv(OutfileName_GEV+'.csv')

                    # Gumble_out_df.to_csv(OutfileName_Gumble+'.csv')

                    pieces = {"GEV": GEV_out_df, "Gumble": Gumble_out_df}
                    result_scenario_model = pd.concat(pieces)

                    Dict_scenario[model] = result_scenario_model
                Dict_scenario_df = pd.concat(Dict_scenario)

                ### From this point on it is added to create the summary of the GEV and gumble distributions
                print("Station and Scenarios are :", station, Scenario)
                print(" 90th percentile mean, 5th, and 95th,  ")
                print(Dict_scenario_df["0.9"].mean(), np.percentile(Dict_scenario_df["0.9"], 5),
                      np.percentile(Dict_scenario_df["0.9"], 95))
                SumArray = [station, Scenario, "0.9", Dict_scenario_df["0.9"].mean(),
                            np.percentile(Dict_scenario_df["0.9"], 5), np.percentile(Dict_scenario_df["0.9"], 95)]
                sumArraySeries.append(SumArray)
                print("Station and Scenarios are :", station, Scenario)
                print(" 99th percentile mean, 5th, and 95th,  ")
                print(Dict_scenario_df["0.99"].mean(), np.percentile(Dict_scenario_df["0.99"], 5),
                      np.percentile(Dict_scenario_df["0.99"], 95))
                SumArray = [station, Scenario, "0.99", Dict_scenario_df["0.99"].mean(),
                            np.percentile(Dict_scenario_df["0.99"], 5), np.percentile(Dict_scenario_df["0.99"], 95)]
                sumArraySeries.append(SumArray)
    gc.collect()
    Summary_df = pd.DataFrame(data=np.array(sumArraySeries),
                              columns=['Station', 'Scenario', 'return_period', 'mean', '5th_percentile',
                                       '95th_percentile'])
print(dict_not_enough_years)
# save the file
Summary_df.to_csv('All_stations_GCM_IDF_12hr_summary.csv')
del Summary_df, sumArraySeries, SumArray, Dict_scenario, Dict_scenario_df, result_scenario_model, (
    file_paths), data, df, precipitation, loc_g, scale_g
gc.collect()


#### PART 3: 6-HR bias corrected GCM projected data

# 6-HR temporally disaggregated bias corrected data
sumArraySeries = []
dict_not_enough_years = {}
all_missing = 0
for si in range(len(merged_StationList)):
    print(si, ":", len(merged_StationList))
    station = merged_StationList[si]
    file_paths = glob(os.path.join(dir_6master_biasCorrected, station + "_bias_corr_6HRmax_df.csv"))

    if len(file_paths) == 0:
        pass
    else:
        # Import the bias corrected rainfall data:
        for file_path in file_paths:
            # file_path = "drive/MyDrive/Hydrology/DAL_CanESM5_2015_SSP245_addition_linear_scaling.csv"
            df_station = pd.read_csv(file_path)
            # Extract the file name from the path
            file_name = os.path.basename(file_path).split('.')[0]

            # Define the start date and the number of days
            start_date = "2015-01-01"
            num_days = 31411
            num_years = 86
            start_year = "2015"

            for Scenario in Scenarios:
                Dict_scenario = {}
                for model in Models:
                    GEV_out_df = pd.DataFrame(columns=['0.9', '0.99'], index=Methods)
                    OutfileName_GEV = station + "_GEV_" + model + "_" + Scenario

                    Gumble_out_df = pd.DataFrame(columns=['0.9', '0.99'], index=Methods)
                    OutfileName_Gumble = station + "_Gumbel_" + model + "_" + Scenario
                    for method in Methods:
                        df_filtered = df_station.loc[(df_station['model'] == model) & \
                                                     (df_station['methods'] == method) & \
                                                     (df_station['scenario'] == Scenario), \
                            ['year', 'Predicted_6hr_pr', "stid"]]

                        # Generate the date range
                        time_index = pd.date_range(start=start_year, freq="YS", periods=num_years)

                        # Create the DataFrame with this index
                        df = pd.DataFrame(index=time_index)
                        df.index.name = "year"
                        df['Predicted_6hr_pr'] = None
                        try:
                            df['Predicted_6hr_pr'] = df_filtered['Predicted_6hr_pr'].values
                        except ValueError:
                            print("Length of values did not match, dataset years < 86")
                        if len(df_filtered['Predicted_6hr_pr'].values) < 86:
                            print("correcting the mismatching length value issue by filling forward")
                            df_filtered['year'] = pd.to_datetime(df_filtered['year'], format='%Y')
                            df_filtered_fullrange = df_filtered.set_index('year').reindex(time_index)
                            # sum the number of missing years
                            n_missing = 86 - len(df_filtered['Predicted_6hr_pr'].values)
                            all_missing = all_missing + n_missing
                            print("Number of missing years measurements", all_missing)
                            # fill in the dictionary
                            not_enough_measures = {
                                all_missing: {'stid': station,
                                              'number_missing_years': n_missing,
                                              'missing_year': df_filtered_fullrange[df_filtered_fullrange[
                                                  'Predicted_6hr_pr'].isnull()].index.year.values,
                                              'scenario': Scenario,
                                              'model': model,
                                              'methods': method}
                            }
                            print(not_enough_measures)
                            dict_not_enough_years.update(not_enough_measures)
                            # fill forward the missing information
                            df_filtered_ff = df_filtered_fullrange.ffill()
                            df['Predicted_6hr_pr'] = df_filtered_ff['Predicted_6hr_pr'].values
                        else:
                            df['Predicted_6hr_pr'] = df_filtered['Predicted_6hr_pr'].values

                        #### this block is for computing the max values
                        data = df
                        precipitation = data.Predicted_6hr_pr

                        shape, loc, scale = gev.fit(precipitation.values, 0)

                        if np.all(precipitation.values == precipitation.values[0]):
                            continue
                        loc_g, scale_g = gumbel_r.fit(precipitation.values)

                        # create vector of years
                        years = np.arange(1.1, 100, 0.1)

                        output_dir = dir
                        Output_filename = station + "_" + model + "_" + method + "_" + Scenario

                        # calcuting the probablities of GEV and Gumble distributions

                        vals = gev.ppf([0.9, 0.99], shape, loc=loc)
                        GEV_out_df.loc[method] = vals

                        vals = gumbel_r.ppf([0.9, 0.99], loc=loc_g)
                        Gumble_out_df.loc[method] = vals

                    # these two lines are comments for memory save

                    # GEV_out_df.to_csv(OutfileName_GEV+'.csv')

                    # Gumble_out_df.to_csv(OutfileName_Gumble+'.csv')

                    pieces = {"GEV": GEV_out_df, "Gumble": Gumble_out_df}
                    result_scenario_model = pd.concat(pieces)

                    Dict_scenario[model] = result_scenario_model
                Dict_scenario_df = pd.concat(Dict_scenario)

                ### From this point on it is added to create the summary of the GEV and gumble distributions
                print("Station and Scenarios are :", station, Scenario)
                print(" 90th percentile mean, 5th, and 95th,  ")
                print(Dict_scenario_df["0.9"].mean(), np.percentile(Dict_scenario_df["0.9"], 5),
                      np.percentile(Dict_scenario_df["0.9"], 95))
                SumArray = [station, Scenario, "0.9", Dict_scenario_df["0.9"].mean(),
                            np.percentile(Dict_scenario_df["0.9"], 5), np.percentile(Dict_scenario_df["0.9"], 95)]
                sumArraySeries.append(SumArray)
                print("Station and Scenarios are :", station, Scenario)
                print(" 99th percentile mean, 5th, and 95th,  ")
                print(Dict_scenario_df["0.99"].mean(), np.percentile(Dict_scenario_df["0.99"], 5),
                      np.percentile(Dict_scenario_df["0.99"], 95))
                SumArray = [station, Scenario, "0.99", Dict_scenario_df["0.99"].mean(),
                            np.percentile(Dict_scenario_df["0.99"], 5), np.percentile(Dict_scenario_df["0.99"], 95)]
                sumArraySeries.append(SumArray)
    gc.collect()
    Summary_df = pd.DataFrame(data=np.array(sumArraySeries),
                              columns=['Station', 'Scenario', 'return_period', 'mean', '5th_percentile',
                                       '95th_percentile'])
print(dict_not_enough_years)
# save the file
Summary_df.to_csv('All_stations_GCM_IDF_6hr_summary.csv')

print("FINISHED")
# Finished

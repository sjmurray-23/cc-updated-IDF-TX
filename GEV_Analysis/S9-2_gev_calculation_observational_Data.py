# -*- coding: utf-8 -*-
"""
This script is writen to compute the GEV and gumble ditribution for  extreme precipitation observational data

"""
## Importing necessary libraries 
import tensorflow as tf
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import tempfile
from scipy import stats
from scipy.stats import genextreme as gev
from scipy.stats import gumbel_r
from glob import glob


###Defining function to compute empirical return lebel
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
dir_master_obsrv = os.path.join(cwd_path, "Data/Station_24HR_Maximums/Observation_24HRmax")

##listing scenarios
Scenarios = ["observation"]
startyears_dic = {}

## reading starting date and the number of years from the stations CSV data
dir_master_ASOS = os.path.join(cwd_path, '/Data/ASOS/')
file_path_ASOS = os.path.join(dir_master_ASOS, 'TX_ASOS_Final.csv')

df_ASOS_stations = pd.read_csv(file_path_ASOS)
merged_StationList = df_ASOS_stations['stid'].tolist()

#######creating the output database and populating it with the projection values 

sumArraySeries = []
for station in merged_StationList:
    file_paths = glob(os.path.join(dir_master_obsrv, station + "_obs_24HRmax_df.csv"))

    for file_path in file_paths:
        df_station = pd.read_csv(file_path)
        # Extract the file name from the path
        file_name = os.path.basename(file_path).split('.')[0]

        Dict_scenario = {}
        out_df = pd.DataFrame(columns=['station', 'method', '0.9', '0.99'])
        OutfileName_GEV = station + "_GEV_obs"

        Gumble_out_df = pd.DataFrame(columns=['0.9', '0.99'])
        OutfileName_Gumble = station + "_Gumbel_obs"

        # Generate the date range
        time_index = pd.date_range(start=df_station.iloc[0]['year'],
                                   periods=1 + df_station.iloc[-1]['year'] - df_station.iloc[0]['year'], freq="Y")

        # Create the DataFrame with this index
        df = pd.DataFrame(index=time_index)
        df.index.name = "year"
        df['pr'] = None
        df['pr'] = df_station['pr'].values
        data = df
        precipitation = data.pr

        shape, loc, scale = gev.fit(precipitation.values, 0)

        if np.all(precipitation.values == precipitation.values[0]):
            continue
        loc_g, scale_g = gumbel_r.fit(precipitation.values)

        vals_gev = gev.ppf([0.9, 0.99], shape, loc=loc)

        vals_gumble = gumbel_r.ppf([0.9, 0.99], loc=loc_g)

        ### From this point on it is added to create the summary of the GEV and gumble distributions

        SumArray = [station, "0.9", vals_gev[0], vals_gumble[0], (vals_gev[0] + vals_gumble[0]) / 2.0]
        sumArraySeries.append(SumArray)
        SumArray = [station, "0.99", vals_gev[1], vals_gumble[1], (vals_gev[1] + vals_gumble[1]) / 2.0]
        sumArraySeries.append(SumArray)

Summary_obs_df = pd.DataFrame(data=np.array(sumArraySeries),
                              columns=['Station', 'return_period', 'GEV', 'Gumble', 'Mean'])

separator = ", "
Summary_obs_df.to_csv('GEV_Analysis/All_stations_obsrv_24hr_summary.csv')

#### from this part foward the code calucaltes the change factors it needs the file name

Summary_GCM_df = pd.read_csv('GEV_Analysis/All_stations_GCM_IDF_24hr_summary.csv')

for station in merged_StationList:
    for RP in ['0.9', '0.99']:
        df_divisors = float(Summary_obs_df.loc[(Summary_obs_df['Station'] == station) & \
                                               (Summary_obs_df['return_period'] == RP), \
            ['Mean']]['Mean'].tolist()[0])

        Summary_GCM_df.loc[(Summary_GCM_df['Station'] == station) & (Summary_GCM_df['return_period'] == float(RP)), \
            ['mean', '5th_percentile', '95th_percentile']] = \
            Summary_GCM_df.loc[(Summary_GCM_df['Station'] == station) & (Summary_GCM_df['return_period'] == float(RP)), \
                ['mean', '5th_percentile', '95th_percentile']] / df_divisors
Summary_GCM_df.to_csv('GEV_Analysis/All_stations_24hr_CF_summary.csv')

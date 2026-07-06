# -*- coding: utf-8 -*-

"""""
AIM: Apply bias correction techniques on LOCA2 climate models at ASOS stations

1 Part code for applying bias correction on climate models at ASOS station data previously prepared 
Extracting GCMs, SSPs, and Observation netcdf daily data at each station to apply bias corrections

Loop setup to open a netcdf file by climate model and then scenario to then access the data to apply bias correction
across all stations, save the bias correction results, then close the netcdfs

Three bias correction techniques: linear scaling, empirical quantile mapping, quantile delta mapping

"""""

import gc
import os
import pickle

import xarray as xr
from cmethods import adjust
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import datetime
import random
from xclim.core import units

# define functions used in script

# function to get the filepaths from a directory
def list_nc_files(directory, mdl):
    return [os.path.join(directory+f"{mdl}/pr/", file) for file in os.listdir(directory+f"{mdl}/pr/") if
            os.path.isfile(os.path.join(directory+f"{mdl}/pr/", file)) and file.endswith('.nc')]


# function to create saving folder pathway for each model and scenario
def create_saving_folders(folder_path):
    if not os.path.exists(folder_path):
        # if here current thread is stopped and the same dir is created in other thread
        # the next line will raise an exception
        os.makedirs(folder_path)
        print("directory made")
    else: print("directory already exists")


# function to test to see if file already exists for station ID, model, ssp, method, and kind
def test_methodFile_exists(saving_path, station_id, mdl, sspnum, method, kind):
    if kind == "+":
        kind = "addition" 
    else:
        kind = "multiplication"
    if method in ["quantile_delta_mapping", "quantile_mapping"]:
        file_path = f'{saving_path}/' + station_id + f'_{mdl}_2015_SSP{sspnum}_nq50{kind}_' + method + '.csv'
    else:
        file_path = f'{saving_path}/' + station_id + f'_{mdl}_2015_SSP{sspnum}_{kind}_' + method + '.csv'
    print(file_path)
    if not os.path.exists(file_path):
        # if here method and kind should continue to perform analysis and save the folder
        return "file does not exist"
    else:
        return "file already exists"


# setup
# initial start time
ist = datetime.datetime.now()
print("Initial script start time: ", ist)

# Set project directory
cwd = os.getcwd()
print(cwd)

# Step 1: load data
# set observation data pathways
OBS_location=cwd+"Data/ASOS/PreparedData/Daily/"
# load ASOS lat lon data
# ASOS station locations path
asos_xy_file = cwd+"Data/ASOS/TX_ASOS_locs_Final.csv" # Final has only the stations that will be used in the analysis i.e., clean data
# Load and prepare ASOS station Lat/Lon
asos_xy=pd.read_csv(asos_xy_file)

# Step 2: prepare for looping across ASOS stations to apply bias correction
# get the station IDs, lat, lons
asos_xy = asos_xy[['stid', 'lat', 'lon']]
# get station id's found in asos_final
sites = asos_xy['stid'].values.tolist()
len(sites)
# set stid as index
asos_xy.set_index('stid', inplace=True)

# set climate scenarios
SSPs = ['245', '585']

# set climate models
models = ['ACCESS-CM2','CanESM5', 'GFDL-ESM4', 'BCC-CSM2-MR', 'EC-Earth3', 'CNRM-CM6-1', 'CNRM-ESM2']

# set variable
variable_OBS = "pr"
variable_mod = "pr"

# set seeds
np.random.seed(10)
random.seed(10)

# loop across cmip models and scenarios
for model in models:
    # Step 3: set folder paths
    # set climate model saving folder path and create
    model_folder = f'BiasCorrection/BiasCorrected_Data/{model}_biascorr'
    create_saving_folders(model_folder)
    for ssp in SSPs:
        print(f"begin extraction of cmip6 model {model} SSP{ssp}")
        # create scenario saving folder
        saving_folder = f'BiasCorrection/BiasCorrected_Data/{model}_biascorr/SSP{ssp}'
        create_saving_folders(saving_folder)
        # create Figure saving folder
        fig_folder = f'BiasCorrection/BiasCorrected_Data/{model}_biascorr/SSP{ssp}/Fig'
        create_saving_folders(fig_folder)
        # list climate simulation historical file
        simh_file = list_nc_files(prj_dir + "Data/LOCA2_CMIP6_6km/historical/", model)
        # list climate simulation projection files
        simp_files = list_nc_files(prj_dir + f"Data/LOCA2_CMIP6_6km/ssp{ssp}/", model)

        # Loop across stations
        for i in range(len(sites)):
            print(i, ":", len(sites), "sites")

            # Step 4: pull station details
            # site id
            stid = sites[i]
            print("station id: ", stid)
            # lat
            stn_lat=asos_xy.loc[stid, 'lat'] # lat goes here
            # long
            stn_lon=asos_xy.loc[stid, 'lon'] # long goes here
            # start date for station
            start = asos_xy.loc[stid, 'new_start_date'] # use new start date post cleaning
            print(stid, ' observation starting date: ', start)

            # Step 5: extract and prepare netcdf data at station
            # start dates cannot go before the simh (GCM historical) start date '1950-01-01',
            # so if it does the start date must set back to default
            if start < '1950-01-01':
                start = '1950-01-01'
            else:
                pass

            # slice data before modeling
            # Observation (ASOS) netcdf
            obsh = xr.open_dataset(OBS_location + "ASOS_precipitation_TXYdims.nc").sel(
                time=slice(start, '2014-12-31')).sel(lat=stn_lat, lon=stn_lon, method='nearest')
            print("completed obsh xarray opening")
            # GCM Historical
            simh = xr.open_dataset(simh_file[0]).sel(lat=stn_lat, lon=(stn_lon+360), method='nearest').sel(
                time=slice(start, '2014-12-31')) # must use [0]
            print("completed simh xarray opening")
            # GCM Projected
            simp = xr.open_mfdataset(simp_files, parallel=True, combine='by_coords', chunks={'time':1825}).sel(
                lat=stn_lat, lon=(stn_lon+360), method='nearest')
            print("completed simp xarray opening")

            # convert units of precip into correct format
            # first convert from rate to absolute values
            obsh_pr = units.rate2amount(obsh.pr)
            simh_pr = units.rate2amount(simh.pr)
            simp_pr = units.rate2amount(simp.pr)
            print("completed pr rate conversion")
            # next convert to in
            obsh_pr = units.convert_units_to(obsh_pr, "in", context="hydro")
            simh_pr = units.convert_units_to(simh_pr, "in", context="hydro")
            simp_pr = units.convert_units_to(simp_pr, "in", context="hydro")
            print("completed pr unit conversion to inches")

            # Step 6: loop through different bias correction methods and apply
            for mth in ["linear_scaling", "quantile_mapping", "quantile_delta_mapping"]:
                for knd in ["+", "*"]:
                    # check if this method already is completed for a station
                    methodFile = test_methodFile_exists(saving_folder, stid, model, ssp, mth, knd)
                    print(methodFile)
                    if methodFile == "file already exists":
                        print("file already exists, skipping method: ", mth, "; kind: ", knd, f"for {stid}")
                    else:
                        print("method: ", mth, "; kind: ", knd)
                        if mth in ["linear_scaling"]:
                            # apply bias correction
                            results = adjust(method=mth,
                                             obs=obsh_pr,
                                             simh=simh_pr,
                                             simp=simp_pr.chunk(dict(time=-1)),
                                             kind=knd,
                                             group="time.month",
                                             max_scaling_factor= 10)

                            print(f"completed bias correction of {stid} for {mth}")

                            # convert to dataframe for saving
                            df = results.to_dataframe()

                            # multiplication kind
                            if knd == "*":
                                # save as a csv
                                df.to_csv(f'{saving_folder}/' + stid + f'_{model}_2015_SSP{ssp}_multiplication_' + mth + '.csv')

                                # plot as figure
                                figure = df.plot().get_figure()
                                figure.savefig(
                                    f'{fig_folder}/' + stid + f'_{model}_2015_SSP{ssp}_multiplication_' + mth + '.png')
                                plt.close()

                                # clean up
                                del results, df, figure
                                gc.collect()
                                print(f"saved results for method: {mth}, kind: {knd}")

                            # addition kind
                            else:
                                # save as a csv
                                df.to_csv(f'{saving_folder}/' + stid + f'_{model}_2015_SSP{ssp}_addition_' + mth + '.csv')

                                # plot as figure
                                figure = df.plot().get_figure()
                                figure.savefig(
                                    f'{fig_folder}/' + stid + f'_{model}_2015_SSP{ssp}_addition_' + mth + '.png')
                                plt.close()

                                # clean up
                                del results, df, figure
                                gc.collect()
                                print(f"saved results for method: {mth}, kind: {knd}")
                        elif mth in ["quantile_mapping", "quantile_delta_mapping"]:
                            # methods have additional parameter: number of quantiles
                            # set number of quantiles
                            nq = 50

                            # apply bias correction
                            results = adjust(method=mth,
                                             obs=obsh_pr,
                                             simh=simh_pr,
                                             simp=simp_pr.chunk(dict(time=-1)),
                                             n_quantiles=nq,
                                             kind=knd)
                            print(f"completed bias correction of {stid} for {mth}")

                            # convert to dataframe for saving
                            df = results.to_dataframe()

                            # multiplication kind
                            if knd == "*":
                                # save as csv
                                df.to_csv(f'{saving_folder}/' + stid + f'_{model}_2015_SSP{ssp}_nq{nq}multiplication_' + mth + '.csv')

                                # plot
                                figure = df.plot().get_figure()
                                figure.savefig(
                                    f'{fig_folder}/' + stid + f'_{model}_2015_SSP{ssp}_nq{nq}multiplication_' + mth + '.png')
                                plt.close()

                                # clean up
                                del results, df, figure
                                gc.collect()
                                print(f"saved results for method: {mth}, kind: {knd}")

                            # addition kind
                            else:
                                # save as csv
                                df.to_csv(f'{saving_folder}/' + stid + f'_{model}_2015_SSP{ssp}_nq{nq}addition_' + mth + '.csv')

                                # plot
                                figure = df.plot().get_figure()
                                figure.savefig(
                                    f'{fig_folder}/' + stid + f'_{model}_2015_SSP{ssp}_nq{nq}addition_' + mth + '.png')
                                plt.close()

                                # clean up
                                del results, df, figure
                                gc.collect()
                                print(f"saved results for method: {mth}, kind: {knd}")
            # release memory
            obsh.close()
            simh.close()
            simp.close()
            del obsh, simh, simp, obsh_pr, simh_pr, simp_pr
            gc.collect()

            print("memory cleared")
            print(f"{model} - SSP{ssp} for station {stid} completed at...")

            # get current time
            ct2 = datetime.datetime.now()
            print(ct2)
            print(i, ":", len(sites), "sites")
    print(f"Bias corrections are complete for {model} all ", len(sites), " sites")
    print("models left: ", len(models) - models.index(model))
    # get current time
    ct3 = datetime.datetime.now()
    print(ct3)

# Finished
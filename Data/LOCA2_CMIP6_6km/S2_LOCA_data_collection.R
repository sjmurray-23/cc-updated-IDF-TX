# Aim: Import LOCA Climate Projection Data from online server

# Data info located at https://loca.ucsd.edu/loca-version-2-for-north-america-ca-jan-2023/
# Data source server located at https://cirrus.ucsd.edu/~pierce/LOCA2/
# CMIP6 downscaled data using LOCA method at 6km res, only electing for regional split - south

# 4 Part script

## Part 1: set local attributes
## Part 2: Extract CMIP6 Historical precipitation ["pr"]
## Part 3: Extract CMIP6 SSP245 precipitation ["pr"]
## Part 4: Extract CMIP6 SSP585 precipitation ["pr"]

############################
# set working directory
wd = getwd()
setwd(wd)
# set saving data path
dir.data = paste0('/Data')
# create new destination folder (i.e., folder)
dest_folder_LOCA = paste0(dir.data, '/LOCA2_CMIP6_6km')
dir.create(paste0(dest_folder_LOCA))

# set url base pathway for regional splits data
loca_url = paste0("https://cirrus.ucsd.edu/~pierce/LOCA2/CONUS_regions_split")

# packages
## base needs
library(lubridate)
library(tidyverse)
## for net cdfs
library(ncdf4)
library(CFtime)

# configure timeout
getOption('timeout')
# rough estimation for sec*(GB*MBconversion)
ts = 1*(2*1000)
ts/60 # number of minutes
options(timeout = ts)
getOption('timeout')

############################

# Download historical data
# Through the GDP
############################
# first set permanent pathway parameters
# experiment
e = "historical"
# ensemble member
"r1i1p1f1"
"r1i1p1f2" # alternative for CNRM-CM6-1, CNRM-ESM2-1, CNRM-CM6-1-HR

# create new directories (i.e., folder)
dir.create(file.path(dest_folder_LOCA, e))

dest_folder_LOCAH = file.path(dest_folder_LOCA, e)

# first set unique pathways
# model
models = c("ACCESS-CM2", "BCC-CSM2-MR", "CanESM5",
           "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3",
           "GFDL-ESM4")
# variable
variables = c("pr")

# create new directories of model names
lapply(models, function(x) if(!dir.exists(file.path(dest_folder_LOCAH, x))) dir.create(file.path(dest_folder_LOCAH, x)))


# track time for download of the entire dataset
start.t <- Sys.time()
start.t
# loop to get the historical data for each variable and identified model
for (v in variables){
  print(v)
  if (v == "pr"){
    for (m in models){
      print(m)
      if (m %in% c("CNRM-CM6-1", "CNRM-ESM2-1", "CNRM-CM6-1-HR")) {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f2"
        # input the unique parameters
        loca_mev_url = paste(loca_url, m, "south/0p0625deg", em, e, v, sep = "/")
        # set file name
        loca_file = paste(v, m, e, em, "1950-2014.LOCA_16thdeg_v20220519.south.nc",sep = ".")
        # combine for final url path
        loca_final_url = paste(loca_mev_url, loca_file, sep = "/")
        # set the destination file pathway and save to directory pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        destfile = file.path(dest_folder, loca_file)
        # download
        print(paste0("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile), download.file(loca_final_url, destfile, mode = "wb"), "File already downloaded")
        print(paste0("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      else {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f1"
        # input the unique parameters
        loca_mev_url = file.path(loca_url, m, "south/0p0625deg", em, e, v)
        # set file name
        loca_file = paste(v, m, e, em, "1950-2014.LOCA_16thdeg_v20220519.south.nc",sep = ".")
        # combine for final url path
        loca_final_url = file.path(loca_mev_url, loca_file)
        # set the destination file pathway and save to directory pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        destfile = file.path(dest_folder, loca_file)
        # download
        print(paste0("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile), download.file(loca_final_url, destfile, mode = "wb"), "File already downloaded")
        print(paste0("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      gc()
    }
  }
  else {
    for (m in models){
      print(m)
      if (m %in% c("CNRM-CM6-1", "CNRM-ESM2-1", "CNRM-CM6-1-HR")) {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f2"
        # input the unique parameters
        loca_mev_url = file.path(loca_url, m, "south/0p0625deg", em, e, v)
        # set file name
        loca_file = paste(v, m, e, em, "1950-2014.LOCA_16thdeg_v20220413.south.nc",sep = ".")
        # combine for final url path
        loca_final_url = file.path(loca_mev_url, loca_file)
        # set the destination file pathway and save to directory pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        destfile = file.path(dest_folder, loca_file)
        # download
        print(paste0("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile), download.file(loca_final_url, destfile, mode = "wb"), "File already downloaded")
        print(paste0("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      else {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f1"
        # input the unique parameters
        loca_mev_url = file.path(loca_url, m, "south/0p0625deg", em, e, v)
        # set file name
        loca_file = paste(v, m, e, em, "1950-2014.LOCA_16thdeg_v20220413.south.nc",sep = ".")
        # combine for final url path
        loca_final_url = file.path(loca_mev_url, loca_file)
        # set the destination file pathway and save to directory pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        destfile = file.path(dest_folder, loca_file)
        # download
        print(paste0("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile), download.file(loca_final_url, destfile, mode = "wb"), "File already downloaded")
        print(paste0("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      gc()
    }
  }
  gc()
}
time_whole_step <- Sys.time() - start.t
time_whole_step

############################

# Download SSP2 RCP4.5 data
# Through the Link
############################
# first set permanent pathway parameters
# experiment
e = "ssp245"
# ensemble member
"r1i1p1f1"
"r1i1p1f2" # alternative for CNRM-CM6-1, CNRM-ESM2-1, CNRM-CM6-1-HR
# create new directories (i.e., folder)
dir.create(file.path(dest_folder_LOCA, e))

dest_folder_LOCAH = file.path(dest_folder_LOCA, e)

# first set unique pathways
# model
models = c("ACCESS-CM2", "BCC-CSM2-MR", "CanESM5",
           "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3",
           "GFDL-ESM4")
# variable
variables = c("pr")
# years
yrs = c("2015-2044", "2045-2074", "2075-2100")
# create new directories of model names
lapply(models, function(x) if(!dir.exists(file.path(dest_folder_LOCAH, x))) dir.create(file.path(dest_folder_LOCAH, x)))


# track time for download of the entire dataset
start.t <- Sys.time()
print(start.t)
# loop to get the historical data for each variable and identified model
for (v in variables){
  print(v)
  ## Precip data pull
  if (v == "pr"){
    for (m in models){
      print(m)
      if (m %in% c("CNRM-CM6-1", "CNRM-ESM2-1", "CNRM-CM6-1-HR")) {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f2"
        # input the unique parameters
        loca_mev_url = paste(loca_url, m, "south/0p0625deg", em, e, v, sep = "/")
        # set file names based on each year (rather than adding another loop that may slow it down)
        loca_file1 = paste(v, m, e, em, yrs[1], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        loca_file2 = paste(v, m, e, em, yrs[2], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        loca_file3 = paste(v, m, e, em, yrs[3], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        # combine for final url pathways
        loca_final_url1 = paste(loca_mev_url, loca_file1, sep = "/")
        loca_final_url2 = paste(loca_mev_url, loca_file2, sep = "/")
        loca_final_url3 = paste(loca_mev_url, loca_file3, sep = "/")
        # set the destination folder pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        # set the file names to be downloaded into all same destination folder
        destfile1 = file.path(dest_folder, loca_file1)
        destfile2 = file.path(dest_folder, loca_file2)
        destfile3 = file.path(dest_folder, loca_file3)
        # download all 3 files into the destination folder
        print(paste("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile1), download.file(loca_final_url1, destfile1, mode = "wb"), paste(loca_file1, "already downloaded"))
        ifelse(!file.exists(destfile2), download.file(loca_final_url2, destfile2, mode = "wb"), paste(loca_file2, "already downloaded"))
        ifelse(!file.exists(destfile3), download.file(loca_final_url3, destfile3, mode = "wb"), paste(loca_file3, "already downloaded"))
        print(paste("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      else {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f1"
        # input the unique parameters
        loca_mev_url = paste(loca_url, m, "south/0p0625deg", em, e, v, sep = "/")
        # set file names based on each year (rather than adding another loop that may slow it down)
        loca_file1 = paste(v, m, e, em, yrs[1], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        loca_file2 = paste(v, m, e, em, yrs[2], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        loca_file3 = paste(v, m, e, em, yrs[3], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        # combine for final url pathways
        loca_final_url1 = paste(loca_mev_url, loca_file1, sep = "/")
        loca_final_url2 = paste(loca_mev_url, loca_file2, sep = "/")
        loca_final_url3 = paste(loca_mev_url, loca_file3, sep = "/")
        # set the destination folder pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        # set the file names to be downloaded into all same destination folder
        destfile1 = file.path(dest_folder, loca_file1)
        destfile2 = file.path(dest_folder, loca_file2)
        destfile3 = file.path(dest_folder, loca_file3)
        # download all 3 files into the destination folder
        print(paste("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile1), download.file(loca_final_url1, destfile1, mode = "wb"), paste(loca_file1, "already downloaded"))
        ifelse(!file.exists(destfile2), download.file(loca_final_url2, destfile2, mode = "wb"), paste(loca_file2, "already downloaded"))
        ifelse(!file.exists(destfile3), download.file(loca_final_url3, destfile3, mode = "wb"), paste(loca_file3, "already downloaded"))
        print(paste("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      gc()
    }
  }
  ## Temp data pull
  else {
    for (m in models){
      print(m)
      if (m %in% c("CNRM-CM6-1", "CNRM-ESM2-1", "CNRM-CM6-1-HR")) {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f2"
        # input the unique parameters
        loca_mev_url = file.path(loca_url, m, "south/0p0625deg", em, e, v)
        # set file names based on each year (rather than adding another loop that may slow it down)
        loca_file1 = paste(v, m, e, em, yrs[1], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        loca_file2 = paste(v, m, e, em, yrs[2], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        loca_file3 = paste(v, m, e, em, yrs[3], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        # combine for final url pathways
        loca_final_url1 = paste(loca_mev_url, loca_file1, sep = "/")
        loca_final_url2 = paste(loca_mev_url, loca_file2, sep = "/")
        loca_final_url3 = paste(loca_mev_url, loca_file3, sep = "/")
        # set the destination folder pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        # set the file names to be downloaded into all same destination folder
        destfile1 = file.path(dest_folder, loca_file1)
        destfile2 = file.path(dest_folder, loca_file2)
        destfile3 = file.path(dest_folder, loca_file3)
        # download all 3 files into the destination folder
        print(paste("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile1), download.file(loca_final_url1, destfile1, mode = "wb"), paste(loca_file1, "already downloaded"))
        ifelse(!file.exists(destfile2), download.file(loca_final_url2, destfile2, mode = "wb"), paste(loca_file2, "already downloaded"))
        ifelse(!file.exists(destfile3), download.file(loca_final_url3, destfile3, mode = "wb"), paste(loca_file3, "already downloaded"))
        print(paste("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      else {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f1"
        # input the unique parameters
        loca_mev_url = file.path(loca_url, m, "south/0p0625deg", em, e, v)
        # set file names based on each year (rather than adding another loop that may slow it down)
        loca_file1 = paste(v, m, e, em, yrs[1], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        loca_file2 = paste(v, m, e, em, yrs[2], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        loca_file3 = paste(v, m, e, em, yrs[3], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        # combine for final url pathways
        loca_final_url1 = paste(loca_mev_url, loca_file1, sep = "/")
        loca_final_url2 = paste(loca_mev_url, loca_file2, sep = "/")
        loca_final_url3 = paste(loca_mev_url, loca_file3, sep = "/")
        # set the destination folder pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        # set the file names to be downloaded into all same destination folder
        destfile1 = file.path(dest_folder, loca_file1)
        destfile2 = file.path(dest_folder, loca_file2)
        destfile3 = file.path(dest_folder, loca_file3)
        # download all 3 files into the destination folder
        print(paste("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile1), download.file(loca_final_url1, destfile1, mode = "wb"), paste(loca_file1, "already downloaded"))
        ifelse(!file.exists(destfile2), download.file(loca_final_url2, destfile2, mode = "wb"), paste(loca_file2, "already downloaded"))
        ifelse(!file.exists(destfile3), download.file(loca_final_url3, destfile3, mode = "wb"), paste(loca_file3, "already downloaded"))
        print(paste("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      gc()
    }
  }
  gc()
}
time_whole_step <- Sys.time() - start.t
print(time_whole_step)

############################

# Download SSP5 RCP8.5 data
# Through the Link
############################
# first set permanent pathway parameters
# experiment
e = "ssp585"
# ensemble member
"r1i1p1f1"
"r1i1p1f2" # alternative for CNRM-CM6-1, CNRM-ESM2-1, CNRM-CM6-1-HR
# create new directories (i.e., folder)
dir.create(file.path(dest_folder_LOCA, e))

dest_folder_LOCAH = file.path(dest_folder_LOCA, e)

# first set unique pathways
# model
models = c("ACCESS-CM2", "BCC-CSM2-MR", "CanESM5",
           "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3",
           "GFDL-ESM4")
# variable
variables = c("pr")
# years
yrs = c("2015-2044", "2045-2074", "2075-2100")
# create new directories of model names
lapply(models, function(x) if(!dir.exists(file.path(dest_folder_LOCAH, x))) dir.create(file.path(dest_folder_LOCAH, x)))


# track time for download of the entire dataset
start.t <- Sys.time()
print(start.t)
# loop to get the historical data for each variable and identified model
for (v in variables){
  print(v)
  ## Precip data pull
  if (v == "pr"){
    for (m in models){
      print(m)
      if (m %in% c("CNRM-CM6-1", "CNRM-ESM2-1", "CNRM-CM6-1-HR")) {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f2"
        # input the unique parameters
        loca_mev_url = paste(loca_url, m, "south/0p0625deg", em, e, v, sep = "/")
        # set file names based on each year (rather than adding another loop that may slow it down)
        loca_file1 = paste(v, m, e, em, yrs[1], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        loca_file2 = paste(v, m, e, em, yrs[2], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        loca_file3 = paste(v, m, e, em, yrs[3], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        # combine for final url pathways
        loca_final_url1 = paste(loca_mev_url, loca_file1, sep = "/")
        loca_final_url2 = paste(loca_mev_url, loca_file2, sep = "/")
        loca_final_url3 = paste(loca_mev_url, loca_file3, sep = "/")
        # set the destination folder pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        # set the file names to be downloaded into all same destination folder
        destfile1 = file.path(dest_folder, loca_file1)
        destfile2 = file.path(dest_folder, loca_file2)
        destfile3 = file.path(dest_folder, loca_file3)
        # download all 3 files into the destination folder
        print(paste("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile1), download.file(loca_final_url1, destfile1, mode = "wb"), paste(loca_file1, "already downloaded"))
        ifelse(!file.exists(destfile2), download.file(loca_final_url2, destfile2, mode = "wb"), paste(loca_file2, "already downloaded"))
        ifelse(!file.exists(destfile3), download.file(loca_final_url3, destfile3, mode = "wb"), paste(loca_file3, "already downloaded"))
        print(paste("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      else {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f1"
        # input the unique parameters
        loca_mev_url = paste(loca_url, m, "south/0p0625deg", em, e, v, sep = "/")
        # set file names based on each year (rather than adding another loop that may slow it down)
        loca_file1 = paste(v, m, e, em, yrs[1], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        loca_file2 = paste(v, m, e, em, yrs[2], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        loca_file3 = paste(v, m, e, em, yrs[3], "LOCA_16thdeg_v20220519.south.nc",sep = ".")
        # combine for final url pathways
        loca_final_url1 = paste(loca_mev_url, loca_file1, sep = "/")
        loca_final_url2 = paste(loca_mev_url, loca_file2, sep = "/")
        loca_final_url3 = paste(loca_mev_url, loca_file3, sep = "/")
        # set the destination folder pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        # set the file names to be downloaded into all same destination folder
        destfile1 = file.path(dest_folder, loca_file1)
        destfile2 = file.path(dest_folder, loca_file2)
        destfile3 = file.path(dest_folder, loca_file3)
        # download all 3 files into the destination folder
        print(paste("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile1), download.file(loca_final_url1, destfile1, mode = "wb"), paste(loca_file1, "already downloaded"))
        ifelse(!file.exists(destfile2), download.file(loca_final_url2, destfile2, mode = "wb"), paste(loca_file2, "already downloaded"))
        ifelse(!file.exists(destfile3), download.file(loca_final_url3, destfile3, mode = "wb"), paste(loca_file3, "already downloaded"))
        print(paste("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      gc()
    }
  }
  ## Temp data pull
  else {
    for (m in models){
      print(m)
      if (m %in% c("CNRM-CM6-1", "CNRM-ESM2-1", "CNRM-CM6-1-HR")) {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f2"
        # input the unique parameters
        loca_mev_url = file.path(loca_url, m, "south/0p0625deg", em, e, v)
        # set file names based on each year (rather than adding another loop that may slow it down)
        loca_file1 = paste(v, m, e, em, yrs[1], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        loca_file2 = paste(v, m, e, em, yrs[2], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        loca_file3 = paste(v, m, e, em, yrs[3], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        # combine for final url pathways
        loca_final_url1 = paste(loca_mev_url, loca_file1, sep = "/")
        loca_final_url2 = paste(loca_mev_url, loca_file2, sep = "/")
        loca_final_url3 = paste(loca_mev_url, loca_file3, sep = "/")
        # set the destination folder pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        # set the file names to be downloaded into all same destination folder
        destfile1 = file.path(dest_folder, loca_file1)
        destfile2 = file.path(dest_folder, loca_file2)
        destfile3 = file.path(dest_folder, loca_file3)
        # download all 3 files into the destination folder
        print(paste("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile1), download.file(loca_final_url1, destfile1, mode = "wb"), paste(loca_file1, "already downloaded"))
        ifelse(!file.exists(destfile2), download.file(loca_final_url2, destfile2, mode = "wb"), paste(loca_file2, "already downloaded"))
        ifelse(!file.exists(destfile3), download.file(loca_final_url3, destfile3, mode = "wb"), paste(loca_file3, "already downloaded"))
        print(paste("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      else {
        # track time for download of the entire dataset
        start.time <- Sys.time()
        print(start.time)
        # these have a different ensemble member id
        em = "r1i1p1f1"
        # input the unique parameters
        loca_mev_url = file.path(loca_url, m, "south/0p0625deg", em, e, v)
        # set file names based on each year (rather than adding another loop that may slow it down)
        loca_file1 = paste(v, m, e, em, yrs[1], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        loca_file2 = paste(v, m, e, em, yrs[2], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        loca_file3 = paste(v, m, e, em, yrs[3], "LOCA_16thdeg_v20220413.south.nc",sep = ".")
        # combine for final url pathways
        loca_final_url1 = paste(loca_mev_url, loca_file1, sep = "/")
        loca_final_url2 = paste(loca_mev_url, loca_file2, sep = "/")
        loca_final_url3 = paste(loca_mev_url, loca_file3, sep = "/")
        # set the destination folder pathway
        dest_folder = file.path(dest_folder_LOCAH, m, v)
        ifelse(!dir.exists(dest_folder), dir.create(dest_folder), "Folder exists already")
        # set the file names to be downloaded into all same destination folder
        destfile1 = file.path(dest_folder, loca_file1)
        destfile2 = file.path(dest_folder, loca_file2)
        destfile3 = file.path(dest_folder, loca_file3)
        # download all 3 files into the destination folder
        print(paste("download starting for model:", m, " variable:", v))
        # check to make sure no file already exists
        ifelse(!file.exists(destfile1), download.file(loca_final_url1, destfile1, mode = "wb"), paste(loca_file1, "already downloaded"))
        ifelse(!file.exists(destfile2), download.file(loca_final_url2, destfile2, mode = "wb"), paste(loca_file2, "already downloaded"))
        ifelse(!file.exists(destfile3), download.file(loca_final_url3, destfile3, mode = "wb"), paste(loca_file3, "already downloaded"))
        print(paste("download ended for model:", m, " variable:", v))
        time_one_step <- Sys.time() - start.time
        print(time_one_step)
      }
      gc()
    }
  }
  gc()
}
time_whole_step <- Sys.time() - start.t
print(time_whole_step)
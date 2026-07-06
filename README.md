# Uncertainty Projection of the Intensity Duration Frequency Parameters Framework Applied for the State of Texas up to 2100
This is the code for the application of “Uncertainty Projection of the Intensity Duration Frequency Parameters Framework Applied for the State of Texas up to 2100.” This study aims to provide a framework that incorporates the changes from precipitation projections and uncertainties into statewide level IDF curves for use by state and local decision-makers using Texas, USA as a case study.

These scripts were developed using Python 3.11 and R 4.3

## Table of Contents

- [Workflow](#Workflow)
- [Function](#Function)
- [Structure](#Structure)
- [Requirements](#Requirements)
- [Acknowledgements](#Acknowledgements)

## Workflow
This project contains several components and is performed in the following order

* **S1: ASOS Data Collection**
* **S2: LOCA Data Collection**
* **S3: ASOS Data Generation**
* **S4: ASOS Data Cleaning and Processing**
* **S5: Bias Correction Application**
* **S6: Prepare Projection Data for Temporal Disaggregation**
* **S7: Temporal Disaggregation**
* **S8: Extreme Value Analysis**
* **S9: Extreme Distribution Fitting**
* **S10: Change Factor and Coefficient of Variance Calculation**
* **S11: Prepare Results for Mapping**
* **S12: Perform Interpolation in ArcGIS**

## Function

| Script                                      | Function                                                                                                                                      |
|---------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|
| S1_ASOS_data_collection.py                  | Download ASOS data from IEM online repository access through API https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?help               |
| S2_LOCA_data_collection.R                   | Download LOCA2 CMIP6 precipitation data from LOCA https://loca.ucsd.edu/loca-version-2-for-north-america-ca-jan-2023/                         |
| S3_Obs_data_generation.py                   | Generate hourly, daily, monthly, yearly, and rain event observational data from ASOS data                                                     |
| S4_Obs_cleaning_prep.py                     | Cleaning and preparing the observational data for bias correction and machine learning temporal disaggregation applications                   |
| S5_Bias_CorrectionTechniques_application.py | Apply array of bias correction techniques to projected LOCA2 precipitation                                                                    |
| S6_BiasCorrected_TemporalDisagg_prep.py     | Prepare bias corrected projected precipitation datasets for temporal disaggregation using machine learning methods                            |
| S7_Temporal_Disaggregation_train_pred.R     | Train and apply temporal disaggregation of 24-hour to 12-hour and 6-hour using random forest                                                  |
| S8_Annual_maxima_generation.py              | Extreme value analysis to calculate annual maximas for 24-hour, 12-hour, and 6-hour precipitation                                             |
| S9-1_gev_calculation_Projection_Data.py     | Compute extreme distribution using Generalized Extreme Value (GEV) and Gumble for extreme precipitation projected datasets (24-, 12-, 6-hour) |
| S9-2_gev_calculation_observational_Data.py  | Calculate the extreme values using Generalized Extreme Value (GEV) and Gumble for extreme precipitation observation datasets (24-hour)        |
| S10_Prepare_results_for_mapping.py          | Prepare the IDFs, CFs, and CVs results for interpolation and mapping                                                                          |
| S11_Calculate_CF_CV.py                      | Calculate the change factors and coefficient of variance on change factors on the 24-hour dataset                                             |


## Structure
The source data is available in the data folder.

The scripts are available under the scripts folder with the corresponding name in the workflow. There is a mixture of Python and R code scripts for different aspects of the analysis, see Requirements section for more details.

The folders for saving and used throughout the project workflow are in the project and structured as laid out in the scripts.

## Requirements
See py-requirements.txt for the dependent python packages.

See r-requirements.R for the dependent R packages.

## Acknowledgements
If using this github repository please site the following publication

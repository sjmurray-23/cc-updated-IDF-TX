# Aim
# Temporal disaggregation of 24-HR precipitation to 12-HR and 6-HR events
################################

# following Lucas 2020 https://doi.org/10.1002/ecm.1422

# 3 PART code:
# Part 1: set local workspace and parameters
# Part 2: build a random forest model to temporally disaggregate precipitation
# Part 3: temporally disaggregate bias corrected GCMs precipitation

################################

# Part 1: build a random forest model to temporally disaggregate precipitation

################################
# Step 1: prepare workspace and directories

# empty memory
rm(list = ls())

# set up current wd 
wd = getwd()
setwd(wd)
# set saving paths 
wd.mdls = c('models/')
wd.fig = c('fig/')
wd.data <- c('model_results/')

# library

## General
library(tidyverse)
library(ggplot2)
library(gridExtra)
## Modelling libraries
library(caret) 
library(lime)
library(pdp)
library(ICEbox)
library(iml)
library(elasticnet)
library(ranger)
library(hrbrthemes)

# Function to get the best tuned model based on model tuning parameters decided when building the model
get_best_result <- function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  return(best_result)
}


# Step 2: create the list of grouped stations previously determined
# based on statewide precipitation patterns

site_groups <- list(
  g1 = c('ELP', 'PEQ', 'GDP', 'FST', 'INK'),
  g2 = c('6R6', 'MAF', 'MRF', 'MDD', 'DUX', 'E38', 'ODO'),
  g3 = c('SNK', 'LRD', 'HRX', 'SJT', 'PYX', 'BPC', 'SWW', 'LBB', 'FTN', 'APY', 
         'HBV', 'PVW', 'DHT', 'BPG', 'OZA', 'GNC', 'DRT', 'MFE', 'AMA'),
  g4 = c('UVA', 'BKS', 'ECU', 'JCT', 'COT', 'ALI', 'HRL', 'EBG', 'CDS', 'HHF', 
         'ABI'),
  g5 = c('BRO', 'COM', 'F05', 'PIL', 'T82', 'PEZ', 'AQO', 'BKD', 'SSF', 'NGP', 
         'BEA', 'CWC', 'NQI', 'NOG', 'HDO', 'SPS', 'CVB', 'BMQ', 'DZB', 'RBO'),
  g6 = c('GRK', 'MKN', 'GOP','EDC', 'ERV', 'MNZ', 'RYW', 'LZZ', 'SEP', 'RPH', 
         'RKP', 'HLR', 'GDJ', 'HYI', 'BAZ', 'CRP', 'MWL', 'SAT', 'ATT', 'AUS'),
  g7 = c('5C1', 'NFW', 'DFW', 'PWG', 'GYB', 'AFW', 'INJ', 'LUD', 'GPM', 
         'RWV', 'FTW', 'JWY', 'DAL', 'DTO', 'XBP', 'GLE', 'FWS', 'GTU', 'RBD', 
         'TPL', '0F2', 'LHB', 'ACT', 'GKY', 'RAS'),
  g8 = c('TRL', 'CLL', '11R', 'HQZ', '3T5', 'ADS', 'GYI', 'VCT', 'JDD', 'PKV', 
         'CRS', 'TKI', 'GVT'),
  g9 = c('GGG', 'SLR', 'PRX', 'JXI', 'TYR', 'PSX', 'PSN', 'JSO', 'GLS', 
         'OSA'),
  g10 = c('AXH', 'UTS', 'MCJ', 'SGR', 'IAH', 'OCH', 'RFI', 'DWH', 'CXO', 'LFK', 
          'BYY'),
  g11 = c('LVJ', 'HOU', '6R3', 'LBX', 'JAS'),
  g12 = c('ORG', 'BPT')
)

# check to make sure there are no miss-typed IDs
setdiff(sites, list_c(site_groups))
setdiff(list_c(site_groups), sites)

################################

# Part 2: build a random forest model to temporally disaggregate precipitation
# combine station locations based on groups
# square root transformation
# train, test split
################################
# Step 1: get & load data filepaths

# Folder holding data
folder = "Temporal disaggregation/"
filepath = paste0(wd, folder)
# list all 6-HR files
files_6hr <- list.files(path = filepath, pattern = "_6_ASOSparameters.csv", 
                        full.names = TRUE)
files_12hr <- list.files(path = filepath, pattern = "_12_ASOSparameters.csv", 
                         full.names = TRUE)

# Create dataframes to hold model metrics
all_model_metrics <- data.frame()

# loop across all site groups to build models
for (g in 1:length(site_groups)){
  print(g)
  
  # Step 2: load files in group
  
  # get vector from list
  group <- site_groups[[g]]
  # create empty vector to put files into
  group_files_6hr <- c()
  group_files_12hr <- c()
  
  # limit to files of interest based on station IDs in group
  for (i in 1:length(group)) {
    # 6-HR files for stid
    f1 = files_6hr[grepl(site_groups[[g]][i], files_6hr)]
    group_files_6hr = c(group_files_6hr, f1)
    # 12-HR files for stid
    f2 = files_12hr[grepl(site_groups[[g]][i], files_12hr)]
    group_files_12hr = c(group_files_12hr, f2)
  }
  rm(i)
  
  # read into one dataframe for 6-HR and for 12-HR
  # 6-HR
  df_6hr <- group_files_6hr %>% lapply(read_csv) %>% 
    bind_rows()
  
  # 12-HR
  # raise exception for group 8 ADS 12-hour does not have any data and causes error
  if (g == 8) {
    # leaving out ADS which is number 6 in list
    df_12hr <- group_files_12hr[c(1:5, 7:13)] %>% 
      lapply(read_csv) %>% bind_rows() 
  } else {
    df_12hr <- group_files_12hr %>% lapply(read_csv) %>% 
      bind_rows() 
  }
  
  # Step 3: select 6-HR dataset parameters for modeling
  print("begin 6-HR events model building")
  
  # select model parameters
  dataset <- df_6hr %>% select(Precipitation_y, DOY, elev,
                                 pr_monthly_T, pr_annual_T, pr_1day_T, 
                                 `pr_+1day_T`, `pr_-1day_T`, `pr_+1monthly_T`, 
                                 `pr_-1monthly_T`) %>% drop_na()
  
  # check the number of observations
  print(paste("6-HR number of samples (without NAs)", nrow(dataset)))
  
  # Step 4: 6-HR Pre process data, train/test split and apply data transformation
  # Set random seed
  set.seed(450)
  
  # make train/test partition
  train_index <- caret::createDataPartition(dataset$Precipitation_y, p = .8,
                                            list = FALSE,
                                            times = 1)
  # train split
  df_train <- dataset[train_index, ]
  # test split
  df_test <- dataset[-train_index, ]
  
  # apply data transformation using square root
  df_train_scaled <- df_train %>% mutate_at(vars(-Precipitation_y), sqrt)
  df_test_scaled <- df_test %>% mutate_at(vars(-Precipitation_y), sqrt)
  
  # get number of features
  n_features <- length(setdiff(names(df_train_scaled), "Precipitation_y"))
  
  ## fold 
  folds <- createFolds(df_train_scaled$Precipitation_y, 
                       k = 5, returnTrain = TRUE)
  trcntrl <- trainControl(index = folds, savePredictions = TRUE, 
                          search = 'random')
  
  # Step 5: Build and train RANDOM FOREST (RANGER) models for 6-HR events
  
  # set model parameter
  list_rf_para <- expand.grid(
    mtry = c(4:n_features) # number of covariates used at each split
    # , from 2 to p by 5
    , splitrule = c('variance') # determines how the decision tree
    # 'variance', 'extratrees'
    , min.node.size = c(3, 5, 10) # min number of data points
    # at a leaf prevent overfitting
    )

  # train model
  mdl_rf <- caret::train(Precipitation_y ~ .
                , data = df_train_scaled
                , method = 'ranger'
                , tuneGrid = list_rf_para
                , trControl = trcntrl
                , na.action = na.omit
                , importance = 'permutation'
                , num.trees = 250
                # , preProcess = "range"
                , metric = "MAE"
                , seed = 916
                )
  
  # Save the model
  saving_file = paste0(wd.mdls, "random_forest_model_group", g, "_6hr.Rds")
  saveRDS(mdl_rf, saving_file)
  
  # Step 6: 6-HR check final model results
  print(mdl_rf$finalModel)
  
  # plot model metric results
  png(filename = paste0(wd.fig, "plot_MAE_mdlgroup", g, "_6hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(mdl_rf, metric = "MAE")
  dev.off()
  png(filename = paste0(wd.fig, "plot_RMSE_mdlgroup", g, "_6hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(mdl_rf, metric = "RMSE")
  dev.off()
  png(filename = paste0(wd.fig, "plot_Rsqr_mdlgroup", g, "_6hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(mdl_rf, metric = "Rsquared")
  dev.off()
  
  # view model predictions of training data
  # extract training data observed values
  obs <- mdl_rf[["trainingData"]][[".outcome"]]
  # extract model predicted values
  pred <- mdl_rf[["finalModel"]][["predictions"]]
  
  # quick plot
  png(filename = paste0(wd.fig, "plot_mdlobsvpred_mdlgroup", g, "_6hr.png"), 
      res = 320,width = 1550, height = 1550)
  plot(obs, pred)
  title(main = paste("Group", g, "Model Obs vs. Pred"))
  dev.off()
  
  # Step 7: 6-HR Check the predictive power on test data 
  test_pred <- predict(mdl_rf, newdata = df_test_scaled)

  # plot
  png(filename = paste0(wd.fig, "plot_Testobsvpred_mdlgroup", g, "_6hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(y = test_pred, x = df_test_scaled$Precipitation_y,
      xlab = "obs", ylab = "pred")
  title(main = paste("Group", g, "Model Obs vs. Pred"))
  dev.off()
  
  # variable importance of each factor, table 2 in Lucas 2020
  imp <- caret::varImp(mdl_rf)
  # plot
  png(filename = paste0(wd.fig, "plot_VarImp_mdlgroup", g, "_6hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(imp, main = paste("Group", g, "Model Variable of Importance"))
  dev.off()
  # save the variable importance table
  save(imp, file = paste0(wd.data,'tbl_model_importance_group', g, '_6hr.Rdata'))
  
  # Calculate the predicted metrics
  test_pred_metrics <- postResample(pred= test_pred, 
                                   obs= df_test_scaled$Precipitation_y)
  df_test_pred_met <- data.frame(Test_RMSE = test_pred_metrics[[1]],
                                Test_Rsquared = test_pred_metrics[[2]],
                                Test_MAE = test_pred_metrics[[3]],
                                Test_SampleSize = length(df_test_scaled$Precipitation_y)
                                , Group = paste0("group", g)
                                , Events = "6hr"
                                )
  
  # get final model results
  best_tune <- get_best_result(mdl_rf)
  best_tune['OOB_Rsqr'] <- mdl_rf$finalModel$r.squared
  best_tune['OOB_MSE'] <- mdl_rf$finalModel$prediction.error
  best_tune['SampleSize'] <- mdl_rf$finalModel$num.samples

  # combine with other dataset
  model_metrics <- cbind(best_tune, df_test_pred_met)
  # add to list of all group output results
  all_model_metrics <- rbind(all_model_metrics, model_metrics)
  
  
  # Step 8: select 12-HR dataset parameters for modeling
  print("begin 12-HR events model building")
  # select model parameters
  dataset <- df_12hr %>% select(Precipitation_y, DOY, elev,
                              pr_monthly_T, pr_annual_T, pr_1day_T, 
                              `pr_+1day_T`, `pr_-1day_T`, `pr_+1monthly_T`, 
                              `pr_-1monthly_T`) %>% drop_na()
  
  # check the number of observation
  print(paste("12-HR number of samples (without NAs)", nrow(dataset)))
  
  # Step 9: 12-HR Pre process data, train/test split and apply data transformation
  
  # Set random seed
  set.seed(450)
  # create train/test partition
  train_index <- caret::createDataPartition(dataset$Precipitation_y, p = .8,
                                           list = FALSE,
                                           times = 1)
  # train split
  df_train <- dataset[train_index, ]
  # test split
  df_test <- dataset[-train_index, ]
  
  # scale data, create the scaling with training data
  df_train_scaled <- df_train %>% mutate_at(vars(-Precipitation_y), sqrt)
  df_test_scaled <- df_test %>% mutate_at(vars(-Precipitation_y), sqrt)
  # get number of features
  n_features <- length(setdiff(names(df_train_scaled), "Precipitation_y"))
  
  ## fold 
  folds <- createFolds(df_train_scaled$Precipitation_y, 
                      k = 5, returnTrain = TRUE)
  trcntrl <- trainControl(index = folds, savePredictions = TRUE, 
                         search = 'random')
  
  # Step 10: Build and train RANDOM FOREST (RANGER) models for 12-HR events

  # set model parameter
  list_rf_para <- expand.grid(
    mtry = c(4:n_features) # number of covariates used at each split
    # , from 2 to p by 5
    , splitrule = c('variance') # determines how the decision tree
    # 'variance', 'extratrees'
    , min.node.size = c(3, 5, 10) # min number of data points
    # at a leaf prevent overfitting
  )
  # train model
  mdl_rf <- caret::train(Precipitation_y ~ .
                         , data = df_train_scaled
                         , method = 'ranger'
                         , tuneGrid = list_rf_para
                         , trControl = trcntrl
                         , na.action = na.omit
                         , importance = 'permutation'
                         , num.trees = 250
                         # , preProcess = "range"
                         , metric = "MAE"
                         , seed = 916
  )
  
  # Save the model
  saving_file = paste0(wd.mdls, "random_forest_model_group", g, "_12hr.Rds")
  saveRDS(mdl_rf, saving_file)
  
  # Step 11: 12-HR check final model results
  print(mdl_rf$finalModel)
  
  # plot model metric results
  png(filename = paste0(wd.fig, "plot_MAE_mdlgroup", g, "_12hr.png"), 
      res = 320, width = 1550, height = 1450)
  plot(mdl_rf, metric = "MAE")
  dev.off()
  png(filename = paste0(wd.fig, "plot_RMSE_mdlgroup", g, "_12hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(mdl_rf, metric = "RMSE")
  dev.off()
  png(filename = paste0(wd.fig, "plot_Rsqr_mdlgroup", g, "_12hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(mdl_rf, metric = "Rsquared")
  dev.off()
  
  # review model predictions
  # extract training data observed values
  obs <- mdl_rf[["trainingData"]][[".outcome"]]
  # extract model predicted values
  pred <- mdl_rf[["finalModel"]][["predictions"]]
  
  # quick plot
  png(filename = paste0(wd.fig, "plot_mdlobsvpred_mdlgroup", g, "_12hr.png"), 
      res = 320,width = 1550, height = 1550)
  plot(obs, pred)
  title(main = paste("Group", g, "Model Obs vs. Pred"))
  dev.off()
  
  # Step 12: 12-HR Check the predictive power on test data 
  test_pred <- predict(mdl_rf, newdata = df_test_scaled)
  # plot
  png(filename = paste0(wd.fig, "plot_Testobsvpred_mdlgroup", g, "_12hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(y = test_pred, x = df_test_scaled$Precipitation_y)
  title(main = paste("Group", g, "Model Obs vs. Pred"),
       xlab = "obs", ylab = "pred")
  dev.off()
  
  # variable importance of each factor, table 2 in Lucas 2020
  imp <- caret::varImp(mdl_rf)
  # plot
  png(filename = paste0(wd.fig, "plot_VarImp_mdlgroup", g, "_12hr.png"), 
      res = 320,width = 1550, height = 1450)
  plot(imp, main = paste("Group", g, "Model Variable of Importance"))
  dev.off()
  
  # save the variable importance table
  save(imp, file = paste0(wd.data,'tbl_model_importance_group', g, '_12hr.Rdata'))
  
  # Calculate the predicted metrics
  test_pred_metrics <- postResample(pred= test_pred, 
                                   obs= df_test_scaled$Precipitation_y)
  df_test_pred_met <- data.frame(Test_RMSE = test_pred_metrics[[1]],
                                Test_Rsquared = test_pred_metrics[[2]],
                                Test_MAE = test_pred_metrics[[3]],
                                Test_SampleSize = length(df_test_scaled$Precipitation_y)
                                , Group = paste0("group", g)
                                , Events = "12hr"
  )
  
  # get final model results
  best_tune <- get_best_result(mdl_rf)
  best_tune['OOB_Rsqr'] <- mdl_rf$finalModel$r.squared
  best_tune['OOB_MSE'] <- mdl_rf$finalModel$prediction.error
  best_tune['SampleSize'] <- mdl_rf$finalModel$num.samples

  # combine with other dataset
  model_metrics <- cbind(best_tune, df_test_pred_met)
  # add to list of all group output results
  all_model_metrics <- rbind(all_model_metrics, model_metrics)
  
  # Step 13: clear memory to prevent errors and slow downs
  # clean environment
  rm(dataset, imp, folds, df_12hr, df_6hr, trcntrl, train_index, 
    mdl_rf, df_test, df_train, df_test_scaled, df_train_scaled, list_rf_para,
    f1, f2, obs, pred, test_pred, group_files_6hr, group_files_12hr)
  gc()
}

# save the model final fit result metrics
write.csv(all_model_metrics, file = paste0(wd.data, 'Final_model_metrics.csv'))

# Loop through the model files to re-create plots that didn't save correctly
model_files <- list.files(path = paste0(wd.mdls), recursive = T,
                          full.names = F)[1:24]
for (mf in model_files) {
  mdl_rf <- readRDS(file = paste0(wd.mdls, mf))
  # create saving name from the model file name
  si <- str_locate(mf, 'group')
  saving_name <- str_sub(mf, si[1], -5)
  # plot model metric results
  png(filename = paste0(wd.fig, "plot_MAE_mdl", saving_name ,".png"),
      res = 320,width = 1550, height = 1450)
  print(plot(mdl_rf, metric = "MAE"))
  dev.off()
  png(filename = paste0(wd.fig, "plot_RMSE_mdl", saving_name ,".png"), 
      res = 320,width = 1550, height = 1450)
  print(plot(mdl_rf, metric = "RMSE"))
  dev.off()
  png(filename = paste0(wd.fig, "plot_Rsqr_mdl", saving_name ,".png"),
      res = 320,width = 1550, height = 1450)
  print(plot(mdl_rf, metric = "Rsquared"))
  dev.off()
  # variable importance of each factor, table 2 in Lucas 2020
  imp <- caret::varImp(mdl_rf)
  # plot
  png(filename = paste0(wd.fig, "plot_VarImp_mdl", saving_name ,".png"),
      res = 320,width = 1300, height = 1550)
  print(plot(imp, main = paste("Group", str_sub(mf, si[2]+1, si[2]+1),
                         "Model Variable of Importance")))
  dev.off()
}


################################

# Part 3: temporally disaggregate bias corrected GCMs daily precipitation

################################
# Step 1: List filepaths of bias corrected GCMs

# Folder holding data
folder = "Temporal disaggregation/"
filepath = paste0(wd, folder)
# list all 6-HR files
bc_files <- list.files(path = filepath, pattern = "MLinput.csv",
                        recursive = TRUE)

# folder for saving prediction figures
wd.fig <- paste0(wd.fig, "BiasCorrect_GCM_results/")

# Loop through the groups to temporally disaggregate Bias Corrected GCMs
for (g in 1:length(site_groups)){
  # check group number, will be used to call model file
  print(g)
  # get vector of station IDs from list
  group <- site_groups[[g]]
  
  # loop through station IDs in group (i.e., model)
  for (i in 1:length(group)) {
    # Step 2: list station bias corrected GCM files
    # get station ID
    stid = group[i]
    print(paste("beginning analysis of", stid))
    # list files for station ID
    stid_files = bc_files[grepl(site_groups[[g]][i], bc_files)]
    
    # loop across filepaths to predict 6-HR and 12-HR precipitation
    for (f in 1:length(stid_files)) {
      # Step 3: load station bias corrected GCM file
      # get file
      bc_file <- stid_files[f]
      # read csv files
      bc_df <- read_csv(paste0(filepath, bc_file))
      
      # Step 4: Final Preparation of the dataframe for model prediction
      
      # drop na values
      bc_df_nona <- bc_df %>% drop_na()
      # select model parameters
      dataset <- bc_df_nona %>% select(DOY, elev, pr_monthly_T, 
                                       pr_annual_T, pr_1day_T, `pr_+1day_T`, 
                                       `pr_-1day_T`, `pr_+1monthly_T`, 
                                       `pr_-1monthly_T`)
      # check the number of observations
      print(paste("Number of samples (without NAs)", nrow(dataset)))
      
      # Step 5: BC GCM Pre processing, train/test split and apply data transformation
      
      # Set random seed
      set.seed(3545)
      
      # scale data
      df_scaled <- dataset %>% mutate_all(sqrt)

      # Step 6: Perform Temporal Disaggregation to 6-HR precipitation
      print("starting 6-HR analysis")
      
      # locate appropriate model
      model_file <- list.files(path = wd.mdls, 
                               pattern = paste0("group", g, "_6hr.Rds"), 
                               full.names = T)
      
      # Load the model
      mdl_rf <- readRDS(file = model_file)
      print(paste("model loaded:", model_file))
      
      # Make the prediction
      bc_pred <- predict(mdl_rf, newdata = df_scaled)
      print("completed 6-HR analysis")
      
      # Step 7: Review the 6-HR prediction results
      print("summay statistics of the prediction results")
      print(bc_pred %>% summary())
      
      # add the predictions to the input dataframe
      bc_df_nona["Predicted_6hr_pr"] <- bc_pred
      # check the difference between the 24-HR and 6-HR
      bc_df_nona <- bc_df_nona %>% mutate(Diff_6hr = pr_1day_T - Predicted_6hr_pr)
      # save the predictions
      write.csv(bc_df_nona,
                file = paste0(wd.data, "predicted_group", g, "_6hr_",
                              file_p, ".csv"))
      print("saved 6-HR predictions")
      
      # plot and save a histogram
      file_p <- str_remove(bc_file, "_MLinput.csv")
      png(filename = paste0(wd.fig, "plot_hist_mdlG", g, "_6hr_",
                             file_p, ".png"), 
           res = 320,width = 1550, height = 1550)
      hist(bc_pred, main = paste("Predicted 6-HR Precipitation Events",
                                 str_replace_all(file_p, pattern = "_", " "),
                                 sep="\n"),
           xlab = "Precipitation (inches)")
      dev.off()

      # check the difference between temporally disaggregated and 24-hr
      # plot to see if any are negative
      png(filename = paste0(wd.fig, "plot_difference_G", g, "_6hr_",
                             file_p, ".png"), 
           res = 320,width = 1550, height = 1550)
      plot(bc_df_nona$Diff_6hr, 
           main = paste("Difference Between 24-HR and Predicted 6-HR Precipitation",
                        str_replace_all(file_p, pattern = "_", " "),
                        sep="\n"),
           sub = "Difference = 24-HR - 6-HR Precipitation",
           ylab = "Precipitation Difference (inches)")
      dev.off()
      
      # plot and save a histogram of GCM predicted vs Observation Predicted
      # extract model predicted
      mdl_pred <- mdl_rf$finalModel$predictions
      # create saving name
      file_p <- str_remove(bc_file, "_MLinput.csv")
      png(filename = paste0(wd.fig, "plot_hist_mdlG", g, "_6hr_obs_",
                             file_p, ".png"), 
           res = 320,width = 1550, height = 1550)
      par(mfrow=c(2,1))
      par(mar=c(0,5,3,3))
      hist(bc_pred, main = paste("Predicted 6-HR Precipitation Events",
                                 str_replace_all(file_p, pattern = "_", " "),
                                 sep="\n"),
           xlab = "",
           ylab = "Frequency for GCM",
           breaks = 15, xlim = c(0, max(bc_pred) + 0.25),
           las =1, xaxt="n",
           col = "deepskyblue3")
      box(col="black")
      par(mar=c(5,5,0,3))
      hist(mdl_pred, main = "",
           xlab = "Precipitation (inches)",
           ylab = "Frequency for Observed",
           breaks = 15, xlim = c(0, max(bc_pred) + 0.25), 
           ylim = c(length(mdl_pred)/3, 0),
           las = 1,
           col = "darkgoldenrod")
      box(col="black")
      dev.off()
      # clear model
      rm(mdl_rf)
      
      # Step 8: Perform Temporal Disaggregation to 12-HR precipitation
      print("starting 12-HR analysis")
      # load appropriate model
      model_file <- list.files(path = wd.mdls, 
                               pattern = paste0("group", g, "_12hr.Rds"), 
                               full.names = T)
      # Load the model
      mdl_rf <- readRDS(file = model_file)
      print(paste("model loaded:", model_file))
      
      # Make the prediction
      bc_pred <- predict(mdl_rf, newdata = df_scaled)
      print("completed 12-HR analysis")
      
      # Step 9: Review the 12-HR prediction results
      print("summay statistics of the prediction results")
      print(bc_pred %>% summary())
      
      # add the predictions to the input dataframe
      bc_df_nona["Predicted_12hr_pr"] <- bc_pred
      # check the difference
      bc_df_nona <- bc_df_nona %>% mutate(Diff_12hr = pr_1day_T - Predicted_12hr_pr)
      # save the predictions
      write.csv(bc_df_nona,
                file = paste0(wd.data, "predicted_Tdisagg_group", g,
                              file_p, ".csv"))
      print("saved temporally disaggregated prediction results")
      
      # plot and save a histogram
      file_p <- str_remove(bc_file, "_MLinput.csv")
      png(filename = paste0(wd.fig, "plot_hist_mdlG", g, "_12hr_",
                             file_p, ".png"), 
           res = 320,width = 1550, height = 1550)
      hist(bc_pred, main = paste("Predicted 12-HR Precipitation Events",
                                 str_replace_all(file_p, pattern = "_", " "),
                                 sep="\n"),
           xlab = "Precipitation (inches)")
      dev.off()
      
      # check the difference between temporally disaggregated and 24-hr
      # plot to see if any are negative
      png(filename = paste0(wd.fig, "plot_difference_G", g, "_12hr_",
                             file_p, ".png"), 
           res = 320,width = 1550, height = 1550)
      plot(bc_df_nona$Diff_12hr, 
           main = paste("Difference Between 24-HR and Predicted 12-HR Precipitation",
                        str_replace_all(file_p, pattern = "_", " "),
                        sep="\n"),
           sub = "Difference = 24-HR - 12-HR Precipitation",
           ylab = "Precipitation Difference (inches)")
      dev.off()
      
      # plot and save a histogram of GCM predicted vs Observation Predicted
      # extract model predicted
      mdl_pred <- mdl_rf$finalModel$predictions
      # create saving name
      file_p <- str_remove(bc_file, "_MLinput.csv")
      png(filename = paste0(wd.fig, "plot_hist_mdlG", g, "_12hr_obs_",
                             file_p, ".png"), 
           res = 320,width = 1550, height = 1550)
      par(mfrow=c(2,1))
      par(mar=c(0,5,3,3))
      hist(bc_pred, main = paste("Predicted 12-HR Precipitation Events",
                                 str_replace_all(file_p, pattern = "_", " "),
                                 sep="\n"),
           xlab = "",
           ylab = "Frequency for GCM",
           breaks = 15, xlim = c(0, max(bc_pred) + 0.25),
           las =1, xaxt="n",
           col = "deepskyblue3")
      box(col="black")
      par(mar=c(5,5,0,3))
      hist(mdl_pred, main = "",
           xlab = "Precipitation (inches)",
           ylab = "Frequency for Observed",
           breaks = 15, xlim = c(0, max(bc_pred) + 0.25), 
           ylim = c(length(mdl_pred)/3, 0),
           las = 1,
           col = "darkgoldenrod")
      box(col="black")
      dev.off()
      # clear model
      rm(mdl_rf)
      print(paste("Completed", f, ":", length(stid_files)))
    }
    print(paste("Completed", i, ":", length(group)))
    gc()
  }
  rm(i)
  gc()
}
################################
# Finished
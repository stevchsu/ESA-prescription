library(reshape2)
library(readxl)
library(data.table)
library(lubridate)
library(zoo)
library(doMC)
library(caret)
library(Metrics)
library(onehot)
library(abind)

# Variable 'look_back_times' determines how many times of records we need to look back on this row (visiting record) of the patient.
look_back_times <- 12

# For LSTM-1, we set Hbc_standard_level to 11. For LSTM-2, we set Hbc_standard_level to 11.2
Hbc_standard_level <- 11.2

core_n <- 10

inverse_rescale_num <- function(scaled_number, numeric_scale) {
  original_number <- scaled_number * diff(numeric_scale) + numeric_scale[1]
  return(original_number)
}

rescale_num <- function(original_number, numeric_scale) {
  scaled_number <- (original_number - numeric_scale[1]) / diff(numeric_scale)
  return(scaled_number)
}

##### Read the dataset and do data-preprocessing #####
{
  # Read dataset
  epo_df <- readRDS("<path to data>")
  
  # Remove columns have NA
  missing_numbers_per_column <- sapply(epo_df, function(X) sum(is.na(X))) 
  remove_columns <- names(missing_numbers_per_column[missing_numbers_per_column > 0])
  epo_df[, remove_columns] <- NULL
  
  # Add new column (PRE_Hbc_to_11): How much Pre_Hbc close to 11
  epo_df$PRE_Hbc_to_11 <- epo_df$PRE_Hbc - Hbc_standard_level
  
  # Add new column (delta_Hbc)
  epo_df$delta_Hbc <- epo_df$POST_Hbc - epo_df$PRE_Hbc
  
  # Sort by patient's ID and first visit date
  epo_df <- epo_df[order(epo_df$NAME, epo_df$PRE_FUNC_DATE), ]
  
  # Difference days of two vistings
  DIF_DAY <- as.integer(epo_df$POST_FUNC_DATE - epo_df$PRE_FUNC_DATE)
  # Add new column (week_n): gap of two vistings in week format

  # calculate week_n
  n = rep(0,nrow(epo_df))
  for (i in 1:nrow(epo_df)) {
    if(DIF_DAY[i]>11&DIF_DAY[i]<=18){
      DIF_DAY[i]=14
    }
    else if(DIF_DAY[i]>18&DIF_DAY[i]<=23){
      DIF_DAY[i]=21
    }
    if((DIF_DAY[i]-1) %% 7!=0){
      n[i] = ((DIF_DAY[i]-1)%/% 7 + 1)
    }
    else{
      n[i] = ((DIF_DAY[i]-1)%/% 7)
    }
  }
  
  epo_df$week_n <- n
  
  # Add new column (EPO_DOSE_per_week)
  epo_df$EPO_DOSE_per_week <- epo_df$EPO_DOSE / epo_df$week_n
  
  # Fill in NA cells from the previous observations
  split_epo_df_by_name <- split(epo_df, epo_df$NAME)
  epo_df <- rbindlist(lapply(split_epo_df_by_name, function(x){
    for(i in 1:ncol(x)){
      if( any(is.na(x[, i])) ){
        x[, i] = na.locf(x[, i], na.rm = F)
        x[, i] = na.locf(x[, i], fromLast = T, na.rm = F)
      }
    }
    return(x)
  }) )
  
  epo_df <- as.data.frame(epo_df)

  # Reorder dataset by 'Patient ID' and 'Visting Date'
  epo_df <- epo_df[order(epo_df$NAME, epo_df$PRE_FUNC_DATE), ]
  
  # Add delta_* column
  first_moment_columns <- c("PRE_Hbc", "Pre_BW", "Post_BW", "dif_Pre_Post_BW", "WBC", "RBC",
                            "Hct", "MCV", "MCH", "MCHC", "Platelet", "TP", "Alb", "GOT", "GPT", "Alk_P", 
                            "Chol", "TG", "LDL_C", "HDL_C", "Sugar_AC", "BUN_B", "BUN_A", "BUN_Next", 
                            "Creatinine", "Uric_acid", "Na", "K", "P", "Ca_2plus", "URR", "Kt_V_G", 
                            "Kt_V_D", "nPCR", "TACurea", "Fe", "UIBC", "Ferritin", "Tran_sat", "Al", 
                            "Mg", "Zn", "i_PTH", "HS_CRP", "EPO_DOSE_per_week")
  
  epo_df <- rbindlist(lapply(split(epo_df, epo_df$NAME1),
                             function(X) {
                               for(col in first_moment_columns) {
                                 X[, paste0("delta_", col)] <- c(0, diff(X[, col]))    # First observations of all patients should be 0
                               }
                               return(X) 
                             }) )
  
  epo_df <- as.data.frame(epo_df)
  
  epo_df <- rbindlist(lapply(split(epo_df, epo_df$NAME1), FUN = function(df_by_name) {
    df_by_name$PRE_Hbc_ori <- df_by_name$PRE_Hbc
    df_by_name$tm1_PRE_Hbc <- sapply(seq(nrow(df_by_name)),
                                     FUN = function(idx) {
                                       if(idx == 1) {
                                         return(NA)
                                         } else {
                                           return(df_by_name[idx-1, "PRE_Hbc"])
                                           }
                                       })
    return(df_by_name)
  }))
  
  ## remove first observation, whose tm1_PRE_Hbc is NA, n=305, remained 22344
  epo_df <- epo_df[!is.na(epo_df$tm1_PRE_Hbc),] 
  
  epo_df <- epo_df[epo_df$week_n == 2 | epo_df$week_n == 3, ] # exclude week_n is not 2 or 3, n=296, remained 22048
  epo_df <- epo_df[-which(epo_df$EPO_DRUG_NO == "2EPO20"),] # exclude EPO_DRUG_NO=2EPO20, n=4495, remained 17553
  
  ##### end of section #####
  
  epo_df <- as.data.frame(epo_df)
  
  # Normalization will be performed later, so I keep the original EPO data for the last output
  original_epo_df <- epo_df
  
  # Directly remove unecessary columns
  epo_df[, c("POST_FUNC_DATE", "POST_Hbc", "EPO_DOSE")] <- NULL
  
  # Setup variable types
  time_column_names <- "PRE_FUNC_DATE"
  categorical_column_names <- names(which(sapply(epo_df, is.factor)))
  group_name <- "NAME1"
  patient_groups <- epo_df[, group_name]
  EPO_DRUG_NO <- epo_df[, "EPO_DRUG_NO"]
  random_effect_categorical_column_names <- "EPO_DRUG_NO"
  fixed_effect_categorical_column_names <- setdiff(categorical_column_names, 
                                                   c(group_name, random_effect_categorical_column_names))
  target_column_name <- "delta_Hbc"
  scalable_columns <- setdiff(colnames(epo_df), c(time_column_names, categorical_column_names, target_column_name, "PRE_Hbc_ori", "tm1_PRE_Hbc"))
  fixed_effect_scalable_column_names <- c("age", "ttc")
  random_effect_scalable_column_names <- setdiff(scalable_columns, 
                                                 fixed_effect_scalable_column_names)
  
  # Get the maximum value of dates
  max_date <- max(epo_df[, time_column_names])
 
  train_idx <- which(epo_df[, time_column_names] <= max_date-months(8))
  valid_idx <- which(epo_df[, time_column_names] <= max_date-months(4) & epo_df[, time_column_names] > max_date-months(8))
  test_idx <- which(epo_df[, time_column_names] <= max_date & epo_df[, time_column_names] > max_date-months(4))
  # Remove date column
  epo_df[, time_column_names] <- NULL
  
  # Use scale of training set as a standard scale
  preprocessParams <- preProcess(epo_df[train_idx, scalable_columns], method=c("range"))
  epo_df[, scalable_columns] <- predict(preprocessParams, epo_df[, scalable_columns])
  
 
  # Do one-hot encoding
  need_one_hot_column_names <- categorical_column_names
  for(need_one_hot_column_name in need_one_hot_column_names) {
    onehot_encoded_column <- predict(onehot(epo_df[, need_one_hot_column_name, drop=F], 
                                            max_levels = max(sapply(epo_df[, need_one_hot_column_name, drop=F], 
                                                                    function(X) length(levels(X)) ))), 
                                     data = epo_df) 
    epo_df <- cbind(epo_df, onehot_encoded_column)
    epo_df[, need_one_hot_column_name] <- NULL
  }
  
  {
    # These code is appending historical data onto the latest record by person
    personal_dat <- split(epo_df, patient_groups)
    personal_dat <- personal_dat[sapply(personal_dat, nrow) > 0]
    
    y_dat <- sapply(X = personal_dat, 
                    FUN = function(df) {
                      df[, target_column_name]
                    })
    y_dat <- unlist(y_dat)
    
    # Request parallel resource
    registerDoMC(core_n)
    
    padding_personal_dat <- foreach(df = personal_dat) %dopar% {
      extracted_end_idx <- seq(nrow(df))
      extracted_start_idx <- extracted_end_idx - look_back_times
      padding_row_numbers <- sapply(X = extracted_start_idx, 
                                    FUN = function(start_idx) ifelse(start_idx <= 0, abs(start_idx) + 1, 0))
      extracted_start_idx[extracted_start_idx <= 0] <- 1
      padding_mtx <- lapply(X = seq(nrow(df)), 
                            FUN = function(idx) {
                              extracted_df <- df[extracted_start_idx[idx]:extracted_end_idx[idx], ]
                              extracted_mtx <- as.matrix(extracted_df)
                              extracted_mtx <- rbind(matrix(0, nrow = padding_row_numbers[idx], ncol = ncol(extracted_mtx)), 
                                                     extracted_mtx)
                              extracted_mtx <- apply(extracted_mtx, 2, as.numeric)
                              return(extracted_mtx)
                            })
      return(padding_mtx)
    }
    
    padding_personal_dat <- unlist(padding_personal_dat, recursive = F)
    padding_personal_dat <- abind::abind(padding_personal_dat, along = 0.5)
    
    # Release parallel resource
    registerDoSEQ()
    
    all_column_names <- colnames(padding_personal_dat[1,,])
    fixed_effect_categorical_idx <- unlist(sapply(fixed_effect_categorical_column_names,
                                                  function(X) grep(paste0("^", X, "=.*$"), all_column_names)))
    fixed_effect_scalable_idx <- unlist(sapply(fixed_effect_scalable_column_names,
                                               function(X) grep(paste0("^", X, "$"), all_column_names)))
    random_effect_categorical_idx <- grep(paste0("^", random_effect_categorical_column_names, "=.*$"), all_column_names) # "EPO_DRUG_NO="
    random_effect_scalable_idx <- unlist(sapply(random_effect_scalable_column_names,
                                                function(X) grep(paste0("^", X, "$"), all_column_names)))
    group_idx <- grep(paste0("^", group_name, "=.*$"), all_column_names)
  }
  
  col_names <- dimnames(padding_personal_dat)[[3]]
  EPO_DOSE_per_week_idx <- which(col_names %in% "EPO_DOSE_per_week")

  # Prepare dataset by their own index (by date)
  ready_dat <- Map(idx = list("train" = c(train_idx), "test" = test_idx, "valid" = valid_idx),
                   f = function(idx) {
                     
                     col_names <- dimnames(padding_personal_dat)[[3]]
                     EPO_DOSE_per_week_idx <- which(col_names %in% "EPO_DOSE_per_week")
                     without_EPO_DOSE_per_week_random_effect_scalable_idx <- setdiff(random_effect_scalable_idx, EPO_DOSE_per_week_idx)
                     
                     return(list(x_random = padding_personal_dat[idx,, c(random_effect_categorical_idx, without_EPO_DOSE_per_week_random_effect_scalable_idx)], 
                                 x_group = padding_personal_dat[idx, look_back_times+1, group_idx], 
                                 x_fixed = padding_personal_dat[idx, look_back_times+1, c(fixed_effect_categorical_idx, fixed_effect_scalable_idx)],
                                 y_delta_hbc = y_dat[idx] ))
                   })
  
  # Get possible dose depends on the enumnation of EPO_DRUG_NO, week_n
  possible_dose <- unique(cbind(EPO_DRUG_NO, epo_df[, c("week_n", "EPO_DOSE_per_week")]))
  possible_dose <- possible_dose[order(possible_dose[, "EPO_DRUG_NO"], possible_dose[, "week_n"], possible_dose[, "EPO_DOSE_per_week"]), ]
  EPO_DOSE_per_week_rescale_range <- preprocessParams$ranges[, "EPO_DOSE_per_week"]
  delta_EPO_DOSE_per_week_rescale_range <- preprocessParams$ranges[, "delta_EPO_DOSE_per_week"]
  
  possible_dose <- lapply(split(possible_dose, possible_dose$EPO_DRUG_NO), 
                          function(X1) {
                            lapply(split(X1, X1$week_n),
                                   function(X2) {
                                     EPO_DOSE_per_week <- X2[, "EPO_DOSE_per_week"]
                                     EPO_DRUG_NO <- as.character(unique(X2[, "EPO_DRUG_NO"]))
                                     original_scale_EPO_DOSE_per_week <- EPO_DOSE_per_week * diff(EPO_DOSE_per_week_rescale_range) + EPO_DOSE_per_week_rescale_range[1]
                                     switch(EPO_DRUG_NO,
                                            "2REC50" = EPO_DOSE_per_week[original_scale_EPO_DOSE_per_week <= 2],
                                            "2NES05" = EPO_DOSE_per_week[original_scale_EPO_DOSE_per_week <= 2],
                                            "2EPO20" = EPO_DOSE_per_week[original_scale_EPO_DOSE_per_week <= 6])
                                   })
                          })
  
  # Simulate EPO dose on every visiting record according their week_n and EPO_DOSE_NO (the type of EPO)
  simulated_dat <- Map(idx = seq(dim(padding_personal_dat)[1]),
                       f = function(idx) {
                         # Get temp data
                         {
                           tmp_x_random <- padding_personal_dat[idx,, c(random_effect_categorical_idx, random_effect_scalable_idx)]
                           tmp_x_group <- padding_personal_dat[idx, look_back_times+1, group_idx]
                           tmp_x_fixed <- padding_personal_dat[idx, look_back_times+1, c(fixed_effect_categorical_idx, fixed_effect_scalable_idx)]
                           tmp_y_delta_hbc <- y_dat[idx]
                         }
                         # Get week_n of this record
                         extracted_week_n <- tmp_x_random[look_back_times+1, "week_n"]
                         # Get EPO_DRUG_NO of this record
                         extracted_EPO_DRUG_NO <- EPO_DRUG_NO[idx]
                         # Get simulated EPO_DOSE_per_week of this record
                         possible_EPO_DOSE_per_week <- possible_dose[[extracted_EPO_DRUG_NO]][[as.character(extracted_week_n)]]
                         # Get replication times
                         replication_times <- length(possible_EPO_DOSE_per_week)
                         # Replicate tmp_x_random
                         tmp_x_random <- replicate(replication_times, tmp_x_random)
                         # Change dimension oders
                         tmp_x_random <- aperm(tmp_x_random, c(3,1,2))
                         # Subsitute EPO_DOSE_per_week with simulated EPO_DOSE_per_week
                         tmp_x_random[, look_back_times+1, "EPO_DOSE_per_week"] <- possible_EPO_DOSE_per_week
                         
                         original_possible_EPO_DOSE_per_week <- inverse_rescale_num(scaled_number = possible_EPO_DOSE_per_week,
                                                                                    numeric_scale = EPO_DOSE_per_week_rescale_range)
                         
                         last_EPO_DOSE_per_week <- unique(tmp_x_random[, look_back_times, "EPO_DOSE_per_week"])
                         original_last_EPO_DOSE_per_week <- inverse_rescale_num(scaled_number = last_EPO_DOSE_per_week, 
                                                                                numeric_scale = EPO_DOSE_per_week_rescale_range)
                         
                         original_possible_delta_EPO_DOSE_per_week <- original_possible_EPO_DOSE_per_week - original_last_EPO_DOSE_per_week
                         scaled_possible_delta_EPO_DOSE_per_week <- rescale_num(original_number = original_possible_delta_EPO_DOSE_per_week, 
                                                                                numeric_scale = delta_EPO_DOSE_per_week_rescale_range)
                         
                         tmp_x_random[, look_back_times+1, "delta_EPO_DOSE_per_week"] <- scaled_possible_delta_EPO_DOSE_per_week
                         
                         # Replicate tmp_x_fixed
                         tmp_x_group <- replicate(replication_times, tmp_x_group)
                         tmp_x_group <- t(tmp_x_group)
                         # Replicate tmp_x_fixed
                         tmp_x_fixed <- replicate(replication_times, tmp_x_fixed)
                         tmp_x_fixed <- t(tmp_x_fixed)
                         # Replicate y_delta_hbc
                         tmp_y_delta_hbc <- replicate(replication_times, tmp_y_delta_hbc)
                         
                         return(list(x_random = tmp_x_random, 
                                     x_group = tmp_x_group, 
                                     x_fixed = tmp_x_fixed,
                                     y_delta_hbc = tmp_y_delta_hbc,
                                     replication_times = replication_times))
                       })
}




##### Setup GPU Enviornment #####
library(reticulate)
use_condaenv(condaenv = 'TF1', required = TRUE)
conda_python(envname = 'TF1', conda = 'auto')

library(keras)
library(tensorflow)
# Set global enviornment variables
Sys.setenv(CUDA_PATH = "/usr/local/cuda")
Sys.setenv(CUDA_ROOT = "/usr/local/cuda")
Sys.setenv(CUDA_HOME="/usr/local/cuda")
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/usr/local/cuda/bin", sep = ":"))
# Set avaiable GPU ID
Sys.setenv(CUDA_VISIBLE_DEVICES = "1")
# Set avaiable GPU memory
gpu_options <- tf$GPUOptions(per_process_gpu_memory_fraction = 0.5, allow_growth = T)
sess <- tf$Session(config = tf$ConfigProto(gpu_options = gpu_options))
# Set Keras backend = Tensorflow
k_set_session(sess)

fixed_data_dim <- dim(ready_dat$train$x_fixed)[2]
random_data_dim <- dim(ready_dat$train$x_random)[3]
group_data_dim <- dim(ready_dat$train$x_group)[2]
paste('fix:', fixed_data_dim, 'random:', look_back_times + 1, random_data_dim, 'group:', group_data_dim, sep=' ')

##### Define and fit model #####
{
  fixed_input <- layer_input(shape = list(fixed_data_dim),
                             dtype = "float32", name = "fixed_data")
  random_input <- layer_input(shape = list(look_back_times + 1, random_data_dim), 
                              dtype = "float32", name = "random_data")
  group_input <- layer_input(shape = list(group_data_dim), 
                             dtype = "float32", name = "group_data")
  # Setup branch
  fixed_branch <- fixed_input %>%
    layer_dense(units = 8, use_bias = T, activation = "tanh",
                input_shape = list(fixed_data_dim)) %>%
    layer_batch_normalization() %>% 
    layer_dense(units = 4, use_bias = T, activation = "linear") %>%
    layer_batch_normalization()
  random_branch <- random_input %>%
    layer_masking(input_shape = list(look_back_times + 1, random_data_dim)) %>%
    layer_lstm(units = 4, return_sequences = F, stateful = F, activation = "tanh") %>% 
    # layer_dropout(rate = 0.1) %>%
    layer_dense(units = 2, use_bias = T, activation = "tanh")
  group_branch <- group_input %>%
    layer_dense(units = 2, use_bias = T, activation = "tanh",
                input_shape = list(group_data_dim)) %>%
    layer_batch_normalization()
  # Concate random_branch and group_branch
  final_random_branch <- layer_concatenate(list(random_branch, group_branch))
  final_random_branch <- final_random_branch %>%
    layer_dense(units = 4, use_bias = T, activation = "tanh") %>%
    layer_batch_normalization()
  # final_output <- layer_add(list(fixed_branch, final_random_branch))
  embedding_branch <- layer_concatenate(list(fixed_branch, final_random_branch))
  final_output <- embedding_branch %>%
    layer_dense(units = 4, use_bias = T, activation = "tanh") %>%
    layer_batch_normalization() %>% 
    layer_dense(units = 1, use_bias = T, activation = "linear")
  model <- keras_model(list(fixed_input, random_input, group_input),
                       final_output)
  model$summary()
  model %>%
    compile(
      loss = 'MSE',
      # optimizer = optimizer_rmsprop(lr = 0.002),
      optimizer = optimizer_adam(lr = 0.01),
      metrics = c('MAE')
    )
  # Fit model
  current_timestamp <- as.integer(Sys.time())
  saved_checkpoints_name <- paste0("checkpoints", current_timestamp, ".h5")
  history <- model %>% fit(x = list(fixed_data = ready_dat$train$x_fixed,
                                     random_data = ready_dat$train$x_random,
                                    group_data = ready_dat$train$x_group),
                           y = ready_dat$train$y_delta_hbc, 
                           validation_split = 0.2,
                           batch_size = 512, epochs = 1000,
                           callbacks = list(
                             callback_model_checkpoint(filepath = saved_checkpoints_name,
                                                       monitor = "val_mean_absolute_error", mode = "min",
                                                       save_best_only = T, verbose = T)), view_metrics = T
  )
}


##### Load fine-tuned model #####
saved_checkpoints_name <- '<name of model object>' 
paste('Loading ckpt', saved_checkpoints_name, sep=' ')
new_model <- load_model_hdf5(saved_checkpoints_name)

# Do prediction on tranining set
pred_result <- new_model %>% predict(list(ready_dat$train$x_fixed, ready_dat$train$x_random, ready_dat$train$x_group))
# Do evaluation on tranining set
round(mae(ready_dat$train$y_delta_hbc, pred_result), 4)
round(mse(ready_dat$train$y_delta_hbc, pred_result), 4)
round(rmse(ready_dat$train$y_delta_hbc, pred_result), 4)

# Do prediction on validation set 
pred_result <- new_model %>% predict(list(ready_dat$valid$x_fixed, ready_dat$valid$x_random, ready_dat$valid$x_group))
# Do evaluation on validation set 
round(mae(ready_dat$valid$y_delta_hbc, pred_result), 4)
round(mse(ready_dat$valid$y_delta_hbc, pred_result), 4)
round(rmse(ready_dat$valid$y_delta_hbc, pred_result), 4)

# Do prediction on test set 
pred_result <- new_model %>% predict(list(ready_dat$test$x_fixed, ready_dat$test$x_random, ready_dat$test$x_group))
# Do evaluation on validation set 
round(mae(ready_dat$test$y_delta_hbc, pred_result), 4)
round(mse(ready_dat$test$y_delta_hbc, pred_result), 4)
round(rmse(ready_dat$test$y_delta_hbc, pred_result), 4)


##### Get simulated data #####
all_x_random <- sapply(simulated_dat, `[[`, "x_random")
simulated_observation_numbers <- sapply(simulated_dat, `[[`, "replication_times")
all_x_random <- abind(all_x_random, along = 1)
all_x_fixed <- sapply(simulated_dat, `[[`, "x_fixed")
all_x_fixed <- do.call("rbind", c(all_x_fixed))
all_x_group <- sapply(simulated_dat, `[[`, "x_group")
all_x_group <- do.call("rbind", c(all_x_group))

# Predict on simulated data
pred_result <- new_model %>% predict(list(all_x_fixed, all_x_random[,,-EPO_DOSE_per_week_idx], all_x_group))
pred_df <- cbind(pred_result, 
                 all_x_random[,look_back_times+1,c("EPO_DOSE_per_week", "week_n")], 
                 rep(seq(simulated_observation_numbers), simulated_observation_numbers),
                 rep(original_epo_df$PRE_Hbc, simulated_observation_numbers) + pred_result)
pred_df <- as.data.frame(pred_df)
colnames(pred_df) <- c("pred_delta_hbc", "EPO_DOSE_per_week", "week_n", "ID", "pred_post_hbc")

# Reverse the scale of EPO_DOSE_per_week and week_n
EPO_DOSE_per_week_rescale_range <- preprocessParams$ranges[, "EPO_DOSE_per_week"]
pred_df$EPO_DOSE_per_week <- inverse_rescale_num(scaled_number = pred_df$EPO_DOSE_per_week, numeric_scale = EPO_DOSE_per_week_rescale_range)
week_n_rescale_range <- preprocessParams$ranges[, "week_n"]
pred_df$week_n <- inverse_rescale_num(scaled_number = pred_df$week_n, numeric_scale = week_n_rescale_range)
dcast_df <- dcast(pred_df, ID~EPO_DOSE_per_week, value.var = "pred_delta_hbc")

# Get the closest EPO_DOSE_per_week depends on different thresholds
threshold_pred_dose <- sapply(seq(11, 11.5, 0.1), function(threshold) {
  closest_to_11_test_d <- lapply(X = split(pred_df, pred_df$ID),
                                 FUN = function(X) {
                                   X[which.min(abs(X$pred_post_hbc - threshold)), ]
                                 })
  closest_to_11_test_d <- data.table::rbindlist(closest_to_11_test_d)
  closest_to_11_test_d <- as.data.frame(closest_to_11_test_d)
  return(closest_to_11_test_d$EPO_DOSE_per_week)
})
colnames(threshold_pred_dose) <- paste0("EPO_DOSE_per_week_cloest_to_", seq(11, 11.5, 0.1))

# Mark the predictions on EPO_DOSE_per_week set not in order
non_order_warning_flags <- sapply(split(pred_df, pred_df$ID),
                                  function(X) {
                                    !(all(order(X$EPO_DOSE_per_week) == order(X$pred_post_hbc)))
                                  })

# Combine original EPO data, non-order flags, predicted EPO_DOSE_per_week on different threshold and corresponding predicting values together
pred_epo_df <- cbind(original_epo_df, "Non-order-prediction" = non_order_warning_flags, threshold_pred_dose, dcast_df)

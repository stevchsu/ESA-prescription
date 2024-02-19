library(Metrics)

# load data
source("R/prepare_data.R")

# load functions
source("R/lib/merf_fit.R")
source("R/lib/merf_predict.R")

# Setup X variables (fixed-effect part)

X_column_name <- c('EPO_DOSE', 'week_n', 'SUCROFER_DOSE') # copied from #26 in AI_model_online/code/experiments/Algorithm/merf/2022.06.14_merf_prequential.r

# Setup Z variables (random-effect part)
Z_column_name <- c('PRE_Hbc', 'tm1_EPO_DOSE', 'tm1_delta_Hbc', 'tm3_delta_Hbc') # copied from #29 in AI_model_online/code/experiments/Algorithm/merf/2022.06.14_merf_prequential.r

# Setup always split variables in ranger
always_split_variables <- c('tm1_delta_Hbc', 'tm2_delta_Hbc', 'tm3_delta_Hbc', 'PRE_Hbc_to_11', 'tm1_PRE_Hbc_to_11', 'tm2_PRE_Hbc_to_11',
                            'tm3_PRE_Hbc_to_11', 'tm1_delta_EPO_DOSE', 'tm2_delta_EPO_DOSE', 'tm3_delta_EPO_DOSE', 'tm1_EPO_DOSE_per_week',
                            'tm2_EPO_DOSE_per_week', 'tm3_EPO_DOSE_per_week', 'tm1_SUCROFER_DOSE', 'tm2_SUCROFER_DOSE', 'tm3_SUCROFER_DOSE')

random_eff = union(Z_column_name, always_split_variables)

# Training MERF
merf_m <- merf_fit(train_d = train_d, 
                   X_column_name = X_column_name, 
                   y_column_name = response, 
                   clusters_column_name = group_column, 
                   Z_column_name = random_eff,
                   max_iterations = 100, 
                   ranger_params = list(always.split.variables = always_split_variables,
                                        num.trees = 1000, seed = 2))

saveRDS(merf_m, "<path of model object>")

# Measuring performance
sapply(X = list("train" = train_d, "test" = test_d, "all" = clean_df),
              FUN = function(dat) {
         pred <- merf_predict(merf_obj = merf_m, test_d = dat,
                              X_column_name = X_column_name, 
                              clusters_column_name = group_column,
                              Z_column_name = random_eff
         )
         mae_result <- round(mae(dat[, response], pred), 4)
         mse_result <- round(mse(dat[,response], pred), 4)
         rmse_result <- round(rmse(dat[, response], pred), 4) 
         return(c(MAE = mae_result, MSE = mse_result, RMSE = rmse_result))
       })

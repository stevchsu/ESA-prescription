library(dplyr)

# load data
source("R/prepare_data_0618.R")

# load functions
source("R/lib/merf_fit.R")
source("R/lib/merf_predict.R")

## Load model
merf_model <- readRDS("<path to model object>")

## Setup X variables (fixed-effect part)
X_column_name <- c('EPO_DOSE', 'week_n', 'SUCROFER_DOSE') # copied from #26 in AI_model_online/code/experiments/Algorithm/merf/2022.06.14_merf_prequential.r

## Setup Z variables (random-effect part)
Z_column_name <- c('PRE_Hbc', 'tm1_EPO_DOSE', 'tm1_delta_Hbc', 'tm3_delta_Hbc') # copied from #29 in AI_model_online/code/experiments/Algorithm/merf/2022.06.14_merf_prequential.r
## Setup always split variables in ranger
always_split_variables <- c('tm1_delta_Hbc', 'tm2_delta_Hbc', 'tm3_delta_Hbc', 'PRE_Hbc_to_11', 'tm1_PRE_Hbc_to_11', 'tm2_PRE_Hbc_to_11',
                            'tm3_PRE_Hbc_to_11', 'tm1_delta_EPO_DOSE', 'tm2_delta_EPO_DOSE', 'tm3_delta_EPO_DOSE', 'tm1_EPO_DOSE_per_week',
                            'tm2_EPO_DOSE_per_week', 'tm3_EPO_DOSE_per_week', 'tm1_SUCROFER_DOSE', 'tm2_SUCROFER_DOSE', 'tm3_SUCROFER_DOSE')

random_eff = union(Z_column_name, always_split_variables)


## add new column "status" to epo_df
epo_df$status <- ""
epo_df$status <- ifelse(epo_df$delta_Hbc >= 0 & epo_df$POST_Hbc < 10.8, "1", epo_df$status)
epo_df$status <- ifelse(epo_df$delta_Hbc >= 0 & epo_df$POST_Hbc >= 10.8 & epo_df$POST_Hbc <= 11.2, "2", epo_df$status)
epo_df$status <- ifelse(epo_df$delta_Hbc >= 0 & epo_df$POST_Hbc > 11.2, "3", epo_df$status)
epo_df$status <- ifelse(epo_df$delta_Hbc < 0 & epo_df$POST_Hbc > 11.2, "4", epo_df$status)
epo_df$status <- ifelse(epo_df$delta_Hbc < 0 & epo_df$POST_Hbc >= 10.8 & epo_df$POST_Hbc <= 11.2, "5", epo_df$status)
epo_df$status <- ifelse(epo_df$delta_Hbc < 0 & epo_df$POST_Hbc < 10.8, "6", epo_df$status)

rbind(n = table(epo_df$status), prob = round(prop.table(table(epo_df$status)), 3))

epo_df$Index <- seq(1:nrow(epo_df))

## select part of columns to check result
tmp <- epo_df %>% dplyr::select(Index, NAME1, PRE_FUNC_DATE, PRE_Hbc, EPO_DOSE, delta_Hbc, POST_FUNC_DATE, POST_Hbc, 
                                status, tm1_PRE_Hbc) 

##### simulate EPO dose #####
epo_df_2weeks <- epo_df[epo_df$week_n == 2,]
epo_df_3weeks <- epo_df[epo_df$week_n == 3,]

simulation_merf <- function(epo_dose, df, model, X_column_name, clusters_column_name, Z_column_name) {
  df$true_EPO_DOSE <- df$EPO_DOSE
  df$EPO_DOSE <- epo_dose
  
  df$pred_delta_Hbc <- merf_predict(merf_obj = model, test_d = df,
                                    X_column_name = X_column_name, 
                                    clusters_column_name = group_column,
                                    Z_column_name = random_eff)
  df$pred_POST_Hbc <- df$PRE_Hbc + df$pred_delta_Hbc
  # df$pred_POST_Hbc_to11 <- df$pred_pred_POST_Hbc - 11
  
  df <- df %>% dplyr::select(Index, CHARTNO, NAME1, week_n, PRE_FUNC_DATE, PRE_Hbc, POST_FUNC_DATE, POST_Hbc, 
                             true_EPO_DOSE, delta_Hbc, status, EPO_DOSE, pred_delta_Hbc, pred_POST_Hbc, tm1_PRE_Hbc)
  return(df)
}

epo_df_3weeks_simulate <- rbind(simulation_merf(epo_dose = 0, df = epo_df_3weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff),
                                simulation_merf(epo_dose = 1, df = epo_df_3weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff), 
                                simulation_merf(epo_dose = 2, df = epo_df_3weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff), 
                                simulation_merf(epo_dose = 3, df = epo_df_3weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff), 
                                simulation_merf(epo_dose = 4, df = epo_df_3weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff), 
                                simulation_merf(epo_dose = 5, df = epo_df_3weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff), 
                                simulation_merf(epo_dose = 6, df = epo_df_3weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff)) %>% 
  arrange(Index)

epo_df_2weeks_simulate <- rbind(simulation_merf(epo_dose = 0, df = epo_df_2weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff),
                                simulation_merf(epo_dose = 1, df = epo_df_2weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff), 
                                simulation_merf(epo_dose = 2, df = epo_df_2weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff), 
                                simulation_merf(epo_dose = 3, df = epo_df_2weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff), 
                                simulation_merf(epo_dose = 4, df = epo_df_2weeks, model = merf_model, 
                                                X_column_name = X_column_name, clusters_column_name = group_column, Z_column_name = random_eff)) %>%
  arrange(Index)
##### end of simulation #####

## select best simulated dose 
epo_3week_best_dose <- epo_df_3weeks_simulate %>% group_by(Index) %>% slice(which.min(abs(pred_POST_Hbc - 11)))
epo_2week_best_dose <- epo_df_2weeks_simulate %>% group_by(Index) %>% slice(which.min(abs(pred_POST_Hbc - 11)))

epo_dose_simulation_result <- rbind(epo_2week_best_dose, epo_3week_best_dose) %>% arrange(Index)


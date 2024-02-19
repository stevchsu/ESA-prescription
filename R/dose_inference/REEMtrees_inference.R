library(dplyr)
library(REEMtree)
library(doMC)

# load data
source("R/prepare_data.R")

## flag of using single REEM tree or Bagging REEM trees
bagging_flag <- T

## Load model
model_path <- ifelse (bagging_flag, yes = "<path to model object>", 
                      no = "<path to model object>")

REEMtrees <- readRDS(model_path) 

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
tmp <- epo_df %>% dplyr::select(Index, NAME1, PRE_FUNC_DATE, PRE_Hbc, EPO_DOSE, delta_Hbc, POST_FUNC_DATE, POST_Hbc, status) 

##### simulate EPO dose #####
epo_df_2weeks <- epo_df[epo_df$week_n == 2,]
epo_df_3weeks <- epo_df[epo_df$week_n == 3,]

simulation_REEMtree <- function(epo_dose, df, model) {
  df$true_EPO_DOSE <- df$EPO_DOSE
  df$EPO_DOSE <- epo_dose
  
  if(bagging_flag) {
    registerDoMC(2)
    all_pred <- foreach(REEMtree = model, .combine = cbind) %dopar% {
      pred <- predict(REEMtree, df, allow.new.levels = TRUE, df[["NAME1"]], EstimateRandomEffects = T)
      return(pred)
    }
    registerDoSEQ()
    df$pred_delta_Hbc <- rowMeans(all_pred)
  } else {
    df$pred_delta_Hbc <- predict(model, df, allow.new.levels = TRUE, df[["NAME1"]], EstimateRandomEffects = T)
  }
  
  df$pred_POST_Hbc <- df$PRE_Hbc + df$pred_delta_Hbc
  
  df <- df %>% dplyr::select(Index, CHARTNO, NAME1, week_n, PRE_FUNC_DATE, PRE_Hbc, POST_FUNC_DATE, POST_Hbc, 
                             true_EPO_DOSE, delta_Hbc, status, EPO_DOSE, pred_delta_Hbc, pred_POST_Hbc, tm1_PRE_Hbc)
  return(df)
}

epo_df_3weeks_simulate <- rbind(simulation_REEMtree(epo_dose = 0, df = epo_df_3weeks, model = REEMtrees),
                                simulation_REEMtree(epo_dose = 1, df = epo_df_3weeks, model = REEMtrees), 
                                simulation_REEMtree(epo_dose = 2, df = epo_df_3weeks, model = REEMtrees), 
                                simulation_REEMtree(epo_dose = 3, df = epo_df_3weeks, model = REEMtrees), 
                                simulation_REEMtree(epo_dose = 4, df = epo_df_3weeks, model = REEMtrees), 
                                simulation_REEMtree(epo_dose = 5, df = epo_df_3weeks, model = REEMtrees), 
                                simulation_REEMtree(epo_dose = 6, df = epo_df_3weeks, model = REEMtrees)) %>% 
  arrange(Index)

epo_df_2weeks_simulate <- rbind(simulation_REEMtree(epo_dose = 0, df = epo_df_2weeks, model = REEMtrees),
                                simulation_REEMtree(epo_dose = 1, df = epo_df_2weeks, model = REEMtrees), 
                                simulation_REEMtree(epo_dose = 2, df = epo_df_2weeks, model = REEMtrees), 
                                simulation_REEMtree(epo_dose = 3, df = epo_df_2weeks, model = REEMtrees), 
                                simulation_REEMtree(epo_dose = 4, df = epo_df_2weeks, model = REEMtrees)) %>%
  arrange(Index)
##### end of simulation #####

## select best simulated dose
epo_2week_best_dose <- epo_df_2weeks_simulate %>% group_by(Index) %>% slice(which.min(abs(pred_POST_Hbc - 11)))
epo_3week_best_dose <- epo_df_3weeks_simulate %>% group_by(Index) %>% slice(which.min(abs(pred_POST_Hbc - 11)))

epo_dose_simulation_result <- rbind(epo_2week_best_dose, epo_3week_best_dose) %>% arrange(Index)

rbind(n = table(epo_dose_simulation_result$status), 
      prob = round(prop.table(table(epo_dose_simulation_result$status)), 4))

library(data.table)
library(zoo)
library(lubridate)
library(doMC)

## Import data
epo_df <- readRDS("<path to data>")

epo_df$SUCRO_DRUG_NO <- NULL

# Sort by patient's ID and visting date
epo_df <- epo_df[order(epo_df$NAME1, epo_df$PRE_FUNC_DATE), ]

# Calculate the number of week which means the duration of EPO visit
DIF_DAY <- as.integer(epo_df$POST_FUNC_DATE - epo_df$PRE_FUNC_DATE)

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
  
# Calculate the dose per week in this EPO visit
epo_df$EPO_DOSE_per_week <- epo_df$EPO_DOSE / epo_df$week_n

# Remove columns have NA
missing_numbers_per_column <- sapply(epo_df, function(X) sum(is.na(X))) 
remove_columns <- names(missing_numbers_per_column[missing_numbers_per_column > 0]) # remove CHARTNO, Dry_Weight, pct_dif_Pre_Post_BW, Ca_P, Vit_D, HbA1c, CTR, CalScore
must_keep_columns <- c("CHARTNO", "EPO_NO")
remove_columns <- setdiff(remove_columns, must_keep_columns)
epo_df[, remove_columns] <- NULL

# Add new column: How much Pre_Hbc close to 11
epo_df$PRE_Hbc_to_11 <- epo_df$PRE_Hbc - 11

# Add delta hbc
epo_df$delta_Hbc <- epo_df$POST_Hbc - epo_df$PRE_Hbc

look_back_times <- 4
core_n <- 8
append_columns <- c("PRE_Hbc", "EPO_DOSE", "EPO_DOSE_per_week", "week_n", 
                    "SUCROFER_DOSE", "delta_Hbc", "PRE_Hbc_to_11","PRE_FUNC_DATE")
NA_indicators <- paste(append_columns, "NA", sep = "_")
preallocated_d <- matrix(NA, nrow = 1, ncol = length(c(append_columns, NA_indicators)) * look_back_times)
colnames(preallocated_d) <- as.character(outer(paste0("tm", seq(look_back_times), "_"), c(append_columns, NA_indicators), FUN = "paste0"))
registerDoMC(core_n)

epo_df <- foreach(patient_record = split(epo_df, epo_df$NAME), .combine = rbind) %dopar% {
  patient_record <- cbind(patient_record, preallocated_d)
  print(unique(patient_record$NAME))
  
  for(row_idx in seq(nrow(patient_record), 1)) {
    fetch_history_idx <- row_idx - 1
    if(fetch_history_idx > 0) {
      
      # There are still something that we can fetch out
      end_idx <- fetch_history_idx
      end_idx <- ifelse(end_idx <= 0, 1, end_idx)
      start_idx <- end_idx - look_back_times + 1
      start_idx <- ifelse(start_idx <= 0, 1, start_idx)
      fetch_record <- patient_record[start_idx:end_idx, append_columns]
    } else {
      fetch_record <- patient_record[0, append_columns]
    }
    
    non_NA_indicator_mtx <- matrix(FALSE, nrow = nrow(fetch_record), ncol = length(NA_indicators))
    colnames(non_NA_indicator_mtx) <- NA_indicators
    fetch_record <- cbind(fetch_record, non_NA_indicator_mtx)
    
    # Add empty data to fill up the row number of fetch_record to look_back_times
    empty_record <- data.frame(matrix(0, nrow=look_back_times-nrow(fetch_record), ncol=length(append_columns)))
    colnames(empty_record) <- append_columns
    NA_indicator_mtx <- matrix(TRUE, nrow = nrow(empty_record), ncol = length(NA_indicators))
    colnames(NA_indicator_mtx) <- NA_indicators
    empty_record <- cbind(empty_record, NA_indicator_mtx)
    
    fetch_record <- rbind(empty_record, fetch_record)
    
    for(i in seq(nrow(fetch_record))) {
      row_i <- fetch_record[i, ]
      colname_names <- paste0("tm", look_back_times-i+1, "_", colnames(row_i))
      patient_record[row_idx, colname_names] <- row_i
    }
  }
  return(patient_record)
}
registerDoSEQ()
epo_df <- as.data.frame(epo_df)


# Append delta EPO_dose_per_week
epo_df$delta_EPO_DOSE <- epo_df$EPO_DOSE - epo_df$tm1_EPO_DOSE
epo_df <- epo_df[!epo_df$tm1_EPO_DOSE_NA,] # exclude first obs of each patient, n=305

append_columns <- c("delta_EPO_DOSE")
NA_indicators <- paste(append_columns, "NA", sep = "_")
preallocated_d <- matrix(NA, nrow = 1, ncol = length(c(append_columns, NA_indicators)) * look_back_times)
colnames(preallocated_d) <- as.character(outer(paste0("tm", seq(look_back_times), "_"), c(append_columns, NA_indicators), FUN = "paste0"))
registerDoMC(core_n)

epo_df$NAME1 <- as.character(epo_df$NAME1)
epo_df$NAME1 <- as.factor(epo_df$NAME1)

epo_df <- foreach(patient_record = split(epo_df, epo_df$NAME), .combine = rbind) %dopar% {
  patient_record <- cbind(patient_record, preallocated_d)
  
  for(row_idx in seq(nrow(patient_record), 1)) {
    fetch_history_idx <- row_idx - 1
    
    if(fetch_history_idx > 0) {
      # There are still something that we can fetch out
      end_idx <- fetch_history_idx
      end_idx <- ifelse(end_idx <= 0, 1, end_idx)
      start_idx <- end_idx - look_back_times + 1
      start_idx <- ifelse(start_idx <= 0, 1, start_idx)
      fetch_record <- patient_record[start_idx:end_idx, append_columns, drop = F]
    } else {
      fetch_record <- patient_record[0, append_columns, drop = F]
    }
    
    non_NA_indicator_mtx <- matrix(FALSE, nrow = nrow(fetch_record), ncol = length(NA_indicators))
    colnames(non_NA_indicator_mtx) <- NA_indicators
    fetch_record <- cbind(fetch_record, non_NA_indicator_mtx)
    
    # Add empty data to fill up the row number of fetch_record to look_back_times
    empty_record <- data.frame(matrix(0, nrow=look_back_times-nrow(fetch_record), ncol=length(append_columns)))
    colnames(empty_record) <- append_columns
    NA_indicator_mtx <- matrix(TRUE, nrow = nrow(empty_record), ncol = length(NA_indicators))
    colnames(NA_indicator_mtx) <- NA_indicators
    empty_record <- cbind(empty_record, NA_indicator_mtx)
    
    fetch_record <- rbind(empty_record, fetch_record)
    
    for(i in seq(nrow(fetch_record))) {
      row_i <- fetch_record[i, ]
      colname_names <- paste0("tm", look_back_times-i+1, "_", colnames(row_i))
      patient_record[row_idx, colname_names] <- row_i
    }
  }
  return(patient_record)
}
registerDoSEQ()
epo_df <- as.data.frame(epo_df)

response <- "delta_Hbc"
group_column <- "NAME1"

not_in_train_columns <- c("EPO_NO", "CHARTNO",
                          "POST_Hbc", "POST_FUNC_DATE") # "PRE_FUNC_DATE",

# Filter the observation what DIF_DAY1 is neither 2 nor 3
epo_df <- epo_df[epo_df$week_n == 2 | epo_df$week_n == 3, ] # exclude week_n is not 2 or 3, n=368, remained 32218; exclude 296 obs, remained 22048
epo_df <- epo_df[-which(epo_df$EPO_DRUG_NO == "2EPO20"),] # exclude EPO_DRUG_NO=2EPO20, n=6238, remained 25979

length(unique(epo_df$NAME1))

clean_df <- epo_df
clean_df[, not_in_train_columns] <- NULL

max_date <- max(epo_df$PRE_FUNC_DATE); 

# testing data 
test_d <- clean_df[epo_df$PRE_FUNC_DATE > max_date - months(6), ]

# training data 
train_d <- clean_df[epo_df$PRE_FUNC_DATE <= max_date - months(6), ]


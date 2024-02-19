library(Metrics)
library(REEMtree)
library(doMC)

# load data
source("R/prepare_data.R")

# Model
numOfTrees <- 80
registerDoMC(20)
t <- Sys.time()
REEMtrees_list <- foreach(i = 1:numOfTrees) %dopar% {
  inBagRows = sample(seq(nrow(train_d)), size = nrow(train_d), replace = T)
  OOBagRows = setdiff(seq(nrow(train_d)), inBagRows)
  
  # Column sampling
  exclusive_columns <- c(group_column)
  predictors <- setdiff(colnames(train_d), c(response, exclusive_columns,"Index","PRE_FUNC_DATE"))
  column_sampling_prop <- 0.5
  sample_predictors <- sample(predictors, 
                              length(predictors) * column_sampling_prop, 
                              replace = F)
  REEMtree_formula <- as.formula(paste(response, paste(sample_predictors, collapse = " + "), sep = " ~ "))
  
  # Note that here we set cp = 0 for higher variances
  baggedREEMTree <- REEMtree(REEMtree_formula,
                             data = train_d[inBagRows, ], random = ~1|NAME1, cv = F,
                             tree.control=rpart.control(cp=0.0001, xval = 1), MaxIterations = 500, verbose = T)
  
  tmp_groups <- baggedREEMTree$EffectModel$groups
  baggedREEMTree$EffectModel <- NULL
  baggedREEMTree$EffectModel <- list()
  baggedREEMTree$EffectModel$groups <- tmp_groups
  baggedREEMTree$data <- NULL
  
  pred <- predict(baggedREEMTree ,
                  newdata = train_d[OOBagRows, ],
                  id = train_d[OOBagRows, group_column], EstimateRandomEffects=T)
  
  OOBagMAE <- mean(abs( train_d[OOBagRows, response] - pred ))
  OOBagRMSE <- rmse(train_d[OOBagRows, response], pred)
  gc()
  return(list(baggedREEMTree = baggedREEMTree, OOBagMAE = OOBagMAE, OOBagRMSE = OOBagRMSE))
}

registerDoSEQ()
print(Sys.time() - t)
REEMtrees <- lapply(REEMtrees_list, `[[`, "baggedREEMTree")

saveRDS(REEMtrees, "<save model object>")

predict(REEMtrees[[1]], newdata = test_d, id = test_d$NAME1, allow.new.levels=T, EstimateRandomEffects = T)

bagging_flag <- F

# RMSE report
sapply(X = list("train" = train_d, "test" = test_d),
       FUN = function(dat) {
         
         # Because the model would use quite lots of memory, we suggest you just use 2 cores here to do predicting
         if(bagging_flag) {
           registerDoMC(2)
           all_pred <- foreach(REEMtree = REEMtrees, .combine = cbind) %dopar% {
             pred <- predict(REEMtree, dat, allow.new.levels = TRUE, dat[["NAME1"]], EstimateRandomEffects = T)
             return(pred)
           }
           registerDoSEQ()
           pred <- rowMeans(all_pred)
         } else {
           pred <- predict(REEMtrees[[1]], newdata = dat, dat[["NAME1"]], allow.new.levels=T, EstimateRandomEffects = T)
         }
         
         mae_result <- round(mae(dat[, response], pred), 4) 
         mse_result <- round(mse(dat[,response], pred), 4)
         rmse_result <- round(rmse(dat[, response], pred), 4) 
         return(c(MAE = mae_result, MSE = mse_result, RMSE = rmse_result))
       })

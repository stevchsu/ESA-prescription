library(MazamaCoreUtils)
library(rlist)
library(ranger)
library(MASS)
library(reticulate)

merf_fit <- function(train_d, 
                     X_column_name, Z_column_name = NULL, y_column_name, clusters_column_name, 
                     max_iterations = 20, ranger_params = list()) {
  

  
  # I cannot find corresponding mathematical function of numpy in R, so I decide to import numpy, lol.
  np <- import("numpy")
  
  X <- train_d[, X_column_name, drop = F]  # fixed-effect
  y <- matrix(train_d[, y_column_name], ncol = 1)
  if(is.null(Z_column_name)) {   # random effect
    Z <- matrix(rep(1, nrow(train_d)), ncol = 1)
  } else {
    Z <- train_d[, Z_column_name, drop = F]
    Z <- as.matrix(Z)
  }
  clusters <- train_d[, clusters_column_name, drop = F]
  clusters <- droplevels(clusters)
  
  stopifnot(nrow(Z) == nrow(X))
  stopifnot(nrow(y) == nrow(X))
  stopifnot(nrow(clusters) == nrow(X))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n_clusters <- dim(unique(clusters))[1]
  n_obs <- length(y)
  q <- dim(Z)[2]  # random effects dimension
  
  # Create a series where cluster_id is the index and n_i is the value
  cluster_counts <- table(clusters)
  
  # Do expensive slicing operations only once
  Z_by_cluster <- list()
  y_by_cluster <- list()
  n_by_cluster <- list()
  I_by_cluster <- list()
  indices_by_cluster <- list()
  
  # TODO: Can these be replaced with groupbys? Groupbys are less understandable than brute force.
  for(cluster_id in names(cluster_counts)) {
    # Find the index for all the samples from this cluster in the large vector
    indices_i <- which(clusters == cluster_id)
    indices_by_cluster[[cluster_id]] <- indices_i
    
    # Slice those samples from Z and y
    Z_by_cluster[[cluster_id]] = as.matrix(Z[indices_i, , drop = F])
    y_by_cluster[[cluster_id]] = y[indices_i, , drop = F]
    
    # Get the counts for each cluster and create the appropriately sized identity matrix for later computations
    n_by_cluster[[cluster_id]] = cluster_counts[cluster_id]
    I_by_cluster[[cluster_id]] = diag(cluster_counts[cluster_id])
  }
  
  # Intialize for EM algorithm
  iteration <- 0
  # Note we are using a dataframe to hold the b_hat because this is easier to index into by cluster_id
  # Before we were using a simple numpy array -- but we were indexing into that wrong because the cluster_ids
  # are not necessarily in order.
  b_hat_df <- matrix(0, nrow = n_clusters, ncol = q)
  rownames(b_hat_df) <- names(cluster_counts)
  sigma2_hat <- 1
  D_hat <- diag(q)
  
  # vectors to hold history
  gll_history <- numeric(0)
  gll_early_stop_threshold <- NA
  b_hat_history <- list()
  b_hat_history <- list.append(b_hat_history, b_hat_df)
  sigma2_hat_history <- list()
  sigma2_hat_history <- list.append(sigma2_hat_history, sigma2_hat)
  D_hat_history <- list()
  D_hat_history <- list.append(D_hat_history, D_hat)
  
  early_stop_flag = FALSE
  
  while(iteration < max_iterations && !early_stop_flag) {
    iteration <- iteration + 1
    logger.debug("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    logger.debug(sprintf("Iteration: %d", iteration))
    logger.debug("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ E-step ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # fill up y_star for all clusters
    y_star = matrix(0, ncol = 1, nrow = nrow(y))
    for(cluster_id in names(cluster_counts)) {
      # Get cached cluster slices
      y_i <- y_by_cluster[[cluster_id]]
      Z_i <- Z_by_cluster[[cluster_id]]
      b_hat_i <- b_hat_df[cluster_id, ]  # used to be ix
      logger.debug(sprintf("E-step, cluster %s, b_hat = %f", cluster_id, b_hat_i))
      indices_i <- indices_by_cluster[[cluster_id]]
      
      # Compute y_star for this cluster and put back in right place
      y_star_i <- y_i - Z_i %*% b_hat_i
      y_star[indices_i, ] <- y_star_i[,1]
    }
      
    
    # check that still one dimensional
    # TODO: Other checks we want to do?
    # assert(length(y_star) == 1
    
    # Do the random forest regression with all the fixed effects features
    dat <- cbind(X, y_star)
    rf <- do.call("ranger", c(list(formula = y_star ~ ., data = dat, keep.inbag = T)))#, ranger_params
    all_pred <- predict(rf, dat, predict.all = T)$predictions
    oob_pred_indicator <- sapply(rf$inbag.counts, `==`, 0)
    oob_pred <- sapply(X = seq(nrow(oob_pred_indicator)), 
                       FUN = function(i) {
                         mean(all_pred[i, oob_pred_indicator[i,]])
                       })
    
    f_hat <- matrix(oob_pred, ncol = 1)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ M-step ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sigma2_hat_sum <- 0
    D_hat_sum <- 0
    
    for(cluster_id in names(cluster_counts)) {
      # Get cached cluster slices
      indices_i <- indices_by_cluster[[cluster_id]]
      y_i <- y_by_cluster[[cluster_id]]
      Z_i <- Z_by_cluster[[cluster_id]]
      n_i <- n_by_cluster[[cluster_id]]
      I_i <- I_by_cluster[[cluster_id]]
      
      # index into f_hat
      f_hat_i <- f_hat[indices_i,,drop=F]
      
      # Compute V_hat_i
      V_hat_i <- Z_i %*% D_hat %*% t(Z_i) + sigma2_hat * I_i
      
      # Compute b_hat_i
      V_hat_inv_i <- np$linalg$pinv(V_hat_i)
      logger.debug(sprintf("M-step, pre-update, cluster %s, b_hat = %f", cluster_id, b_hat_df[cluster_id, ]))
      b_hat_i <- D_hat %*% t(Z_i) %*% V_hat_inv_i %*% (y_i - f_hat_i)
      logger.debug(sprintf("M-step, post-update, cluster %s, b_hat = %f", cluster_id, b_hat_df[cluster_id, ]))
      
      # Compute the total error for this cluster
      eps_hat_i <- y_i - f_hat_i - (Z_i %*% b_hat_i)
      
      logger.debug("------------------------------------------")
      logger.debug(sprintf("M-step, cluster %s", cluster_id))
      logger.debug(sprintf("error squared for cluster = %f", t(eps_hat_i) %*% eps_hat_i))
      
      # Store b_hat for cluster both in numpy array and in dataframe
      # Note this HAS to be assigned with loc, otw whole df get erroneously assigned and things go to hell
      b_hat_df[cluster_id, ] <- b_hat_i
      logger.debug(sprintf(
        "M-step, post-update, recalled from db, cluster %s, b_hat = %f", cluster_id, b_hat_df[cluster_id,])
      )
      
      # Update the sums for sigma2_hat and D_hat. We will update after the entire loop over clusters
      sigma2_hat_sum <- sigma2_hat_sum + as.numeric(t(eps_hat_i) %*% eps_hat_i + 
                                                      sigma2_hat * (n_i - sigma2_hat * sum(diag(V_hat_inv_i)))
                                                    )
      D_hat_sum <- D_hat_sum + np$outer(b_hat_i, b_hat_i) + (
        D_hat - D_hat %*% t(Z_i) %*% V_hat_inv_i %*% Z_i %*% D_hat
      ) 
    }
    
    # Normalize the sums to get sigma2_hat and D_hat
    sigma2_hat <- (1.0 / n_obs) * sigma2_hat_sum
    D_hat <- (1.0 / n_clusters) * D_hat_sum
    
    logger.debug(sprintf("b_hat = %f", b_hat_df))
    logger.debug(sprintf("sigma2_hat = %f", sigma2_hat))
    logger.debug(sprintf("D_hat = %f", D_hat))
    
    # Store off history so that we can see the evolution of the EM algorithm
    b_hat_history <- list.append(b_hat_history, b_hat_df)
    sigma2_hat_history <- list.append(sigma2_hat_history, sigma2_hat)
    D_hat_history <- list.append(D_hat_history, D_hat)
    
    # Generalized Log Likelihood computation to check convergence
    gll <- 0
    for(cluster_id in names(cluster_counts)) {
      # Get cached cluster slices
      indices_i <- indices_by_cluster[[cluster_id]]
      y_i = y_by_cluster[[cluster_id]]
      Z_i = Z_by_cluster[[cluster_id]]
      I_i = I_by_cluster[[cluster_id]]
      
      # Slice f_hat and get b_hat
      f_hat_i <- f_hat[indices_i, ]
      R_hat_i <- sigma2_hat * I_i
      b_hat_i <- b_hat_df[cluster_id, ]
      
      # Numerically stable way of computing log(det(A))
      logdet_D_hat = np$linalg$slogdet(D_hat)[[2]]
      logdet_R_hat_i = np$linalg$slogdet(R_hat_i)[[2]]
      
      gll <- gll + (
        as.numeric(t(y_i - f_hat_i - Z_i %*% b_hat_i) %*% np$linalg$pinv(R_hat_i) %*% (y_i - f_hat_i - Z_i %*% b_hat_i)) + 
          as.numeric(t(b_hat_i) %*% np$linalg$pinv(D_hat) %*% b_hat_i) + 
          logdet_D_hat + 
          logdet_R_hat_i
      )
    }
    
    message(sprintf("GLL is %f at iteration %d.", gll, iteration))
    gll_history <- c(gll_history, gll)
    
    # Early Stopping. This code is entered only if the early stop threshold is specified and
    # if the gll_history array is longer than 1 element, e.g. we are past the first iteration.
    if(!is.na(gll_early_stop_threshold) && (length(gll_history) > 1)) {
      curr_threshold <- np.abs((gll - gll_history[length(gll_history)-1]) / gll_history[length(gll_history)-1])
      logger.debug(sprintf("stop threshold = %f", curr_threshold))
      
      if(curr_threshold < gll_early_stop_threshold) {
        message(sprintf("Gll %f less than threshold %f, stopping early ...", gll, curr_threshold))
        early_stop_flag = TRUE
      }
    }
  }
  
  merf_obj <- list(rf = rf, trained_b = b_hat_df, cluster_counts = cluster_counts)
  class(merf_obj) <- "MERF object"
  return(merf_obj)
}



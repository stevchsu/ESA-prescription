library(ranger)


merf_predict <- function(merf_obj, test_d, X_column_name, Z_column_name = NULL, clusters_column_name) {
  
  if(!("rf" %in% names(merf_obj))) {
    stop("This merf instance is not fitted yet. Call 'merf_fit' with appropriate arguments before using this function")
  }
  
  X <- test_d[, X_column_name, drop = F]
  # y <- matrix(test_d[, y_column_name], ncol = 1)
  if(is.null(Z_column_name)) {
    Z <- matrix(rep(1, nrow(test_d)), ncol = 1)
  } else {
    Z <- test_d[, Z_column_name, drop = F]
    Z <- as.matrix(Z)
  }
  clusters <- test_d[, clusters_column_name, drop = F]
  clusters <- droplevels(clusters)
  
  # Apply random forest to all
  y_hat <- predict(merf_obj$rf, X)$predictions
  
  # Apply random effects correction to all known clusters. Note that then, by default, the new clusters get no
  # random effects correction -- which is the desired behavior.
  for(cluster_id in names(merf_obj$cluster_counts)) {
    indices_i <- which(clusters == cluster_id)
    
    # If cluster doesn't exist in test data that's ok. Just move on.
    if(length(indices_i) == 0) {
      next
    }
      
    # If cluster does exist, apply the correction.
    b_i <- merf_obj$trained_b[cluster_id, ]
    Z_i <- Z[indices_i, , drop = F]
    y_hat[indices_i] <- y_hat[indices_i] + Z_i %*% b_i
  }
  
  return(y_hat)

}

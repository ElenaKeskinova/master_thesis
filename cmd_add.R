# function to add points to cmdscale, with Nystrom approximation:
cmd_add = function(D_train, d_new, V ,lambda){
  # - lambda: vector of eigenvalues from MDS
  # - V: matrix of eigenvectors (columns are eigenvectors)
  # - D_train: distance matrix among training points (n x n)
  # - d_new: distances from new point to all training points (vector of length n)
  
  # Step 1: Compute centering terms
  n <- nrow(D_train)
  D2_train <- D_train^2
  mean_rows <- rowMeans(D2_train)
  grand_mean <- mean(D2_train)
  
  # Step 2: Center the new point's distances
  d2_new <- d_new^2
  mean_new <- rowMeans(d2_new)
  K_new <- -0.5 * (d2_new - mean_rows - mean_new + grand_mean)
  
  # Step 3: Project using Nyström formula
  
  y_new <- t(t(K_new %*% (V)) / lambda)
  return(y_new)
}

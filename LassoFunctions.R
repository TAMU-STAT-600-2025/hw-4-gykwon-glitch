# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y) {
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  # [ToDo] Center and scale X
  n <- nrow(X)
  p <- ncol(X)
  Xmeans <- colMeans(X)
  X_centered <- X - matrix(Xmeans, n, p, byrow = TRUE)

  weights <- sqrt(colSums(X_centered^2) / n)
  Xtilde <- sweep(X_centered, 2, weights, "/")
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda) {
  S <- sign(a) * max(abs(a) - lambda, 0)
  return(S)
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda) {
  n <- nrow(Xtilde)
  return((1 / (2 * n)) * sum((Ytilde - Xtilde %*% beta)^2) + lambda * sum(abs(beta)))
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001) {
  # [ToDo]  Check that n is the same between Xtilde and Ytilde
  if ((nrow(Xtilde) != length(Ytilde))) {
    stop("Check for compatibility of dimensions between Xtilde and Ytilde")
  }
  # [ToDo]  Check that lambda is non-negative
  if (lambda < 0) {
    stop("Check that lambda is non-negative")
  }
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  # [ToDo]  Check for starting point beta_start.
  # If none supplied, initialize with a vector of zeros.
  if (is.null(beta_start)) {
    beta <- numeric(p)
  }
  # If supplied, check for compatibility with Xtilde in terms of p
  else {
    if (length(beta_start) != p) stop("check for compatibility with Xtilde")
    beta <- as.numeric(beta_start)
  }
  # [ToDo]  Coordinate-descent implementation.
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  fobj_current <- lasso(Xtilde, Ytilde, beta, lambda)
  diff <- Inf
  while (diff > eps) {
    # r : residual
    r <- as.numeric(Ytilde - (Xtilde %*% beta))

    for (j in 1:p) {
      xj <- Xtilde[, j]
      # a : vector to take soft
      # partial residual: r + xj * beta[j]
      a <- as.numeric(beta[j] + (crossprod(Xtilde[ , j], r) / n))
      beta_j_new <- soft(a, lambda)

      # residual update so that r = Y - Xbeta stays consistent
      r <- r + xj * (beta[j] - beta_j_new)
      beta[j] <- beta_j_new
    }
    # calculate new objective function
    fobj_new <- lasso(Xtilde, Ytilde, beta, lambda)
    # difference between objective functions
    diff <- abs(fobj_current - fobj_new)
    # update
    fobj_current <- fobj_new
  }
  # Return
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  fmin <- fobj_current
  return(list(beta = beta, fmin = fmin))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001) {
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if ((nrow(Xtilde) != length(Ytilde))) {
    stop("Check for compatibility of dimensions between Xtilde and Ytilde")
  }
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  if (!is.null(lambda_seq)){
    lambda_seq <- lambda_seq[which(lambda_seq >= 0)]
    if (length(lambda_seq) == 0 ){
      warning("none of the supplied values satisfy the requirement. proceed as if the values were not supplied")
      lambda_max <- max(abs(1 / n * crossprod(Xtilde, Ytilde)))
      lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
    }
  }
  else {
    # If lambda_seq is not supplied, calculate lambda_max
    lambda_max <- max(abs(1 / n * crossprod(Xtilde, Ytilde)))
    # (the minimal value of lambda that gives zero solution),
    # and create a sequence of length n_lambda as
    lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  # sort 
  lambda_seq <- sort(lambda_seq, decreasing = TRUE)
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda
  # (make sure supplied eps is carried over).
  # Use warm starts strategy discussed in class for setting the starting values.
  lambda_len <- length(lambda_seq)
  beta_mat <- matrix(0, p, lambda_len)
  fmin_vec <- rep(0, lambda_len)
  
  beta_fmin_list <- fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[1], beta_start = beta_mat[ , 1], eps = eps)
  beta_mat[ , 1] <- beta_fmin_list$beta
  fmin_vec[1] <- beta_fmin_list$fmin
  # loop from 2 if lambda_len >= 2
  if(lambda_len >= 2) {
    for(i in 2:lambda_len) {
      # warm start : start from beta_mat[ , i - 1]
      beta_fmin_list <- fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], beta_start = beta_mat[ , (i - 1)], eps = eps)
      beta_mat[ , i] <- beta_fmin_list$beta
      fmin_vec[i] <- beta_fmin_list$fmin
    }
  }
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X, Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001) {
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  standardize <- standardizeXY(X, Y)
  Xtilde <- standardize$Xtilde
  Ytilde <- standardize$Ytilde 
  Ymean <- standardize$Ymean 
  Xmeans <- standardize$Xmeans
  weights <- standardize$weights
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  fLs <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  beta_mat <- sweep(fLs$beta_mat, 1, weights, "/") # beta_mat[j] = beta[j](p*1) / weights(p*1)
  beta0_vec <- Ymean - as.vector(crossprod(Xmeans, beta_mat))
  lambda_seq <- fLs$lambda_seq
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X, Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001) {
  
  n <- nrow(X)
  # [ToDo] Fit Lasso on original data using fitLASSO
  fit <- fitLASSO(X, Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  lambda_seq <- fit$lambda_seq
  L <- length(lambda_seq)
  beta_mat  <- fit$beta_mat
  beta0_vec <- fit$beta0_vec
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if(is.null(fold_ids)){
    if (k <= 1) stop("k must be >= 2 when fold_ids is NULL.")
    # split equally + random order
    order <- sample.int(n)
    fold_ids <- ((seq_len(n) - 1) %% k) + 1
    fold_ids <- fold_ids[order]
  }
  else {
    k <- max(fold_ids)
  }
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  cv_fold <- matrix(0, nrow = k, ncol = L)
  for (fold in seq_len(k)) {
    idx_tr <- fold_ids != fold # train index
    idx_va <- !idx_tr # validation index
    
    Xtr <- X[idx_tr, , drop = FALSE]
    Ytr <- Y[idx_tr]
    Xva <- X[idx_va, , drop = FALSE]
    Yva <- Y[idx_va]
    
    # refitting with the same lambda_seq
    fit_f <- fitLASSO(Xtr, Ytr, lambda_seq = lambda_seq, n_lambda = L, eps = eps)
    num_val <- length(Yva)
    
    # predict
    B0   <- matrix(rep(fit_f$beta0_vec, each = num_val),
                   nrow = num_val, ncol = L, byrow = FALSE)
    Yhat <- B0 + Xva %*% fit_f$beta_mat
    Ymat <- matrix(Yva, nrow = num_val, ncol = L)
    
    cv_fold[fold, ] <- colMeans((Ymat - Yhat)^2)
  }
  
  # cv summary statistics
  cvm  <- colMeans(cv_fold)
  cvse <- apply(cv_fold, 2, stats::sd) / sqrt(k)
  
  # [ToDo] Find lambda_min
  # lambda_min (min-CV rule)
  idx_min    <- which.min(cvm)
  lambda_min <- lambda_seq[idx_min]
  
  # [ToDo] Find lambda_1SE
  # lambda_1se (1-SE rule): Biggest lambda satisfying CV <= minCV + SE(min)
  crit       <- cvm[idx_min] + cvse[idx_min]
  idx_1se    <- which(cvm <= crit)[1]   # lambda_seq is sorted from max to min
  lambda_1se <- lambda_seq[idx_1se]

  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

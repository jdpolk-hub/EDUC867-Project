#Function 1: Fit multiple linear regression using matrices
#Inputs: y = outcome vector; X = predictor matrix without intercept
#Outputs: coefficients, SEs, residual variance, fitted values, t, p, F, R2
#Assumptions: numeric complete data, full-rank matrix, n > p

fit_matrix_reg <- function(y, X) {
  
  #Check inputs
  if (!is.numeric(y) || !is.numeric(as.matrix(X))) {
    stop("Outcome and predictors must be numeric.")
  }
  
  X <- as.matrix(X)
  
  if (anyNA(y) || anyNA(X) ||
      any(!is.finite(y)) || any(!is.finite(X))) {
    stop("Missing or non-finite values are not allowed.")
  }
  
  n <- length(y)
  
  if (nrow(X) != n) {
    stop("Outcome and predictors must have the same number of observations.")
  }
  
  #Add intercept
  X <- cbind(Intercept = 1, X)
  p <- ncol(X)
  
  if (n <= p) {
    stop("Sample size must exceed the number of coefficients.")
  }
  
  if (qr(X)$rank < p) {
    stop("Predictor matrix is rank-deficient.")
  }
  
  #Matrix-based coefficients
  XtX_inv <- solve(t(X) %*% X)
  coefficients <- XtX_inv %*% t(X) %*% y
  
  #Fitted values and residuals
  fitted_values <- X %*% coefficients
  residuals <- y - fitted_values
  
  #Residual variance and standard errors
  SSE <- sum(residuals^2)
  residual_variance <- SSE / (n - p)
  standard_errors <- sqrt(diag(residual_variance * XtX_inv))
  
  #Inferential statistics
  t_statistics <- as.vector(coefficients) / standard_errors
  p_values <- 2 * pt(-abs(t_statistics), df = n - p)
  
  SST <- sum((y - mean(y))^2)
  
  if (SST <= 0) {
    stop("Outcome must have nonzero variance.")
  }
  
  R_squared <- 1 - SSE / SST
  F_statistic <- ((SST - SSE) / (p - 1)) / (SSE / (n - p))
  
  #Return results
  list(
    coefficients = coefficients,
    standard_errors = standard_errors,
    residual_variance = residual_variance,
    fitted_values = fitted_values,
    t_statistics = t_statistics,
    p_values = p_values,
    F_statistic = as.numeric(F_statistic),
    R_squared = as.numeric(R_squared)
  )
}


#################################################################################

#Function 2: Predict outcomes using a fitted matrix regression
#Inputs: model = fitted result; new_X = new predictor matrix
#Output: predicted outcomes
#Assumptions: numeric complete data, same predictor order as fitting

predict_matrix_reg <- function(model, new_X) {
  
  new_X <- as.matrix(new_X)
  
  if (!is.numeric(new_X)) {
    stop("New predictors must be numeric.")
  }
  
  if (anyNA(new_X) || any(!is.finite(new_X))) {
    stop("Missing or non-finite values are not allowed.")
  }
  
  if (ncol(new_X) != nrow(model$coefficients) - 1) {
    stop("Incorrect number of predictor columns.")
  }
  
  #Add intercept
  new_X <- cbind(Intercept = 1, new_X)
  
  #Calculate predictions
  predictions <- new_X %*% model$coefficients
  
  return(as.vector(predictions))
}


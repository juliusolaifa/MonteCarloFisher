mc_fisher_cluster <- function(
    X,
    Z,
    beta,
    family,
    dispersion,
    nsim
) {
  
  # Mean at Gamma = 0
  eta <- drop(X %*% beta)
  mu  <- family$linkinv(eta)
  
  # First simulation establishes score dimension and names
  y <- family$simulate(
    mu         = mu,
    dispersion = dispersion
  )
  
  U0 <- score_glmm(
    y          = y,
    X          = X,
    Z          = Z,
    beta       = beta,
    family     = family,
    dispersion = dispersion
  )
  
  U <- matrix(
    NA_real_,
    nrow = nsim,
    ncol = length(U0),
    dimnames = list(NULL, names(U0))
  )
  
  U[1, ] <- U0
  
  # Remaining simulations
  if (nsim > 1L) {
    
    for (r in 2:nsim) {
      
      y <- family$simulate(
        mu         = mu,
        dispersion = dispersion
      )
      
      U[r, ] <- score_glmm(
        y          = y,
        X          = X,
        Z          = Z,
        beta       = beta,
        family     = family,
        dispersion = dispersion
      )
    }
  }
  
  # Monte Carlo mean score
  score_mean <- colMeans(U)
  
  # MC standard error of the mean score
  score_mcse <- apply(U, 2, sd) / sqrt(nsim)
  
  # Center scores
  U_centered <- sweep(
    U,
    MARGIN = 2,
    STATS  = score_mean,
    FUN    = "-"
  )
  
  # Information for this cluster
  information <- crossprod(U_centered) / nsim
  
  list(
    information = information,
    score_mean  = score_mean,
    score_mcse  = score_mcse
  )
}

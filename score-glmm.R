score_glmm <- function(
    y,
    X,
    Z,
    beta,
    family,
    dispersion = NULL
) {
  
  # Linear predictor at Gamma = 0
  eta <- drop(X %*% beta)
  
  # Mean through the inverse link
  mu <- family$linkinv(eta)
  
  # Family/link-specific derivatives
  deriv <- family$derivatives(
    y          = y,
    mu         = mu,
    dispersion = dispersion
  )
  
  # General GLMM score components
  comp <- subject_score_components(
    X_i = X,
    Z_i = Z,
    g_i = deriv$g,
    h_i = deriv$h
  )
  
  # Covariance-parameter score
  U_gamma <- gamma_score(comp$A)
  
  # Assemble beta and Gamma scores
  U <- c(
    comp$U_beta,
    U_gamma
  )
  
  # Add dispersion score when the family has one
  if (!is.null(deriv$u_dispersion)) {
    
    U_dispersion <- sum(deriv$u_dispersion)
    
    U <- c(
      U,
      setNames(
        U_dispersion,
        family$dispersion_name
      )
    )
  }
  
  U
}

effective_fim <- function(
    object,
    parameters = NULL
) {
  
  if (!inherits(object, "glmm_fim")) {
    stop("object must inherit from class 'glmm_fim'.")
  }
  
  I <- object$I_full
  
  # By default, covariance parameters are the
  # parameters of interest
  if (is.null(parameters)) {
    
    interest <- object$gamma_idx
    
  } else {
    
    if (!all(parameters %in% object$parameters)) {
      stop("Unknown parameter name supplied.")
    }
    
    interest <- match(
      parameters,
      object$parameters
    )
  }
  
  # Everything else is nuisance
  nuisance <- setdiff(
    seq_along(object$parameters),
    interest
  )
  
  # Information for parameters of interest
  I_gg <- I[
    interest,
    interest,
    drop = FALSE
  ]
  
  # No nuisance parameters
  if (length(nuisance) == 0L) {
    return(I_gg)
  }
  
  # Nuisance information
  I_ll <- I[
    nuisance,
    nuisance,
    drop = FALSE
  ]
  
  # Cross-information
  I_gl <- I[
    interest,
    nuisance,
    drop = FALSE
  ]
  
  I_lg <- I[
    nuisance,
    interest,
    drop = FALSE
  ]
  
  # Efficient information:
  #
  # I_{gamma.lambda}
  #   = I_gg - I_gl I_ll^{-1} I_lg
  #
  I_eff <-
    I_gg -
    I_gl %*%
    solve(I_ll, I_lg)
  
  I_eff
}

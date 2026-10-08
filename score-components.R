subject_score_components <- function(X_i, Z_i, g_i, h_i) {
  
  # Fixed-effect score:
  # U_beta,i = X_i^T g_i
  U_beta <- drop(crossprod(X_i, g_i))
  
  if (!is.null(colnames(X_i))) {
    names(U_beta) <- colnames(X_i)
  } else {
    names(U_beta) <- paste0("beta", seq_along(U_beta))
  }
  
  # Random-effect score:
  # S_i = Z_i^T g_i
  S <- drop(crossprod(Z_i, g_i))
  
  # Random-effect Hessian:
  # H_i = Z_i^T diag(h_i) Z_i
  H <- crossprod(Z_i, h_i * Z_i)
  
  # Boundary covariance-score building block:
  # A_i = H_i + S_i S_i^T
  A <- H + tcrossprod(S)
  
  list(
    U_beta = U_beta,
    S      = S,
    H      = H,
    A      = A
  )
}

gamma_score <- function(A) {
  
  q <- nrow(A)
  
  if (ncol(A) != q) {
    stop("A must be square.")
  }
  
  n_gamma <- q * (q + 1) / 2
  
  U_gamma <- numeric(n_gamma)
  nm      <- character(n_gamma)
  
  k <- 1L
  
  for (j in seq_len(q)) {
    for (i in seq_len(j)) {
      
      if (i == j) {
        U_gamma[k] <- 0.5 * A[i, j]
      } else {
        U_gamma[k] <- A[i, j]
      }
      
      nm[k] <- paste0("g", i, j)
      
      k <- k + 1L
    }
  }
  
  names(U_gamma) <- nm
  
  U_gamma
}

# ============================================================
# Negative binomial: vectorized score calculation
# ============================================================

nbinom_score_batch <- function(Y, X, Z, mu, theta) {
  stopifnot(is.matrix(Y), is.matrix(X), is.matrix(Z),
            ncol(Y) == nrow(X), nrow(X) == nrow(Z),
            length(mu) == ncol(Y), length(theta) == 1L,
            is.finite(theta), theta > 0)
  B <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(X)
  q <- ncol(Z)

  # sweep is explicit about recycling across observations (columns).
  G <- sweep(Y, 2L, mu, "-")
  G <- sweep(G, 2L, theta / (theta + mu), "*")

  H <- sweep(Y, 2L, theta, "+")
  H <- sweep(H, 2L, -theta * mu / (theta + mu)^2, "*")

  Ub <- G %*% X
  S  <- G %*% Z

  n_gamma <- q * (q + 1L) / 2L
  Ug <- matrix(0, nrow = B, ncol = n_gamma)
  gamma_names <- character(n_gamma)
  k <- 0L
  for (b in seq_len(q)) {
    for (a in seq_len(b)) {
      k <- k + 1L
      # One matrix-vector product for all B Hessian contributions.
      h_ab <- drop(H %*% (Z[, a] * Z[, b]))
      Ug[, k] <- (h_ab + S[, a] * S[, b]) * if (a == b) 0.5 else 1
      gamma_names[k] <- paste0("g", a, b)
    }
  }

  # theta score; construct a B x m matrix without ambiguous recycling.
  D <- digamma(Y + theta) - digamma(theta) + log(theta) + 1
  D <- sweep(D, 2L, log(theta + mu), "-")
  # Subtract (Y + theta) / (theta + mu), column by column.
  D <- D - sweep(Y + theta, 2L, theta + mu, "/")
  Ut <- rowSums(D)

  beta_names <- colnames(X)
  if (is.null(beta_names)) beta_names <- paste0("beta", seq_len(p))
  colnames(Ub) <- beta_names
  colnames(Ug) <- gamma_names
  U <- cbind(Ub, Ug, theta = Ut)
  storage.mode(U) <- "double"
  U
}


# ============================================================
# Negative binomial family specification
# ============================================================

nbinom_log <- function() {
  
  list(
    family = "Negative binomial",
    link   = "log",
    
    # mu = g^{-1}(eta)
    linkinv = function(eta) {
      exp(eta)
    },
    
    # Simulate one cluster conditional on Gamma = 0
    simulate = function(mu, dispersion) {
      rnbinom(
        n    = length(mu),
        size = dispersion,
        mu   = mu
      )
    },
    
    # Observation-level derivatives
    derivatives = function(y, mu, dispersion) {
      
      theta <- dispersion
      
      # d log f / d eta
      g <- theta * (y - mu) / (theta + mu)
      
      # d^2 log f / d eta^2
      h <- -theta * mu * (theta + y) /
        (theta + mu)^2
      
      # d log f / d theta
      u_dispersion <-
        digamma(y + theta) -
        digamma(theta) +
        log(theta) + 1 -
        log(theta + mu) -
        (y + theta) / (theta + mu)
      
      list(
        g            = g,
        h            = h,
        u_dispersion = u_dispersion
      )
    },

    # Vectorized Monte Carlo score calculation
    score_batch = nbinom_score_batch,
    
    dispersion      = TRUE,
    dispersion_name = "theta"
  )
}

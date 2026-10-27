# ============================================================
# Compound Poisson-Gamma: vectorized score calculation
# ============================================================

cp_score_batch <- function(
    Y, X, Z, mu, phi, power,
    dispersion_score = NULL
) {

  stopifnot(
    is.matrix(Y),
    is.matrix(X),
    is.matrix(Z),
    is.numeric(Y),
    ncol(Y) == nrow(X),
    nrow(X) == nrow(Z),
    length(mu) == ncol(Y),
    all(is.finite(mu)), all(mu > 0),
    length(phi) == 1L,
    is.finite(phi), phi > 0,
    length(power) == 1L,
    is.finite(power),
    power > 1, power < 2
  )

  B <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)

  # First derivative with respect to eta
  G <- sweep(Y, 2L, mu, "-")
  G <- sweep(
    G, 2L, mu^(1 - power) / phi, "*"
  )

  # Second derivative with respect to eta
  H <- sweep(
    Y, 2L, (1 - power) * mu^(1 - power), "*"
  )

  H <- sweep(
    H, 2L, (2 - power) * mu^(2 - power), "-"
  )

  H <- H / phi

  # Fixed-effect score
  Ub <- G %*% X

  # Random-effect intermediate score
  S <- G %*% Z

  # Covariance parameter scores at Gamma = 0
  n_gamma <- q * (q + 1L) / 2L

  Ug <- matrix(0, B, n_gamma)
  gamma_names <- character(n_gamma)

  k <- 0L

  for (b in seq_len(q)) {
    for (a in seq_len(b)) {

      k <- k + 1L

      h_ab <- drop(
        H %*% (Z[, a] * Z[, b])
      )

      Ug[, k] <- (
        h_ab + S[, a] * S[, b]
      ) * if (a == b) 0.5 else 1

      gamma_names[k] <- paste0("g", a, b)
    }
  }

  beta_names <- colnames(X)

  if (is.null(beta_names)) {
    beta_names <- paste0("beta", seq_len(p))
  }

  colnames(Ub) <- beta_names
  colnames(Ug) <- gamma_names

  # Dispersion score requires a separate implementation
  if (is.null(dispersion_score)) {
    U <- cbind(Ub, Ug)
  } else {
    D <- dispersion_score(Y, mu, phi, power)

    stopifnot(
      is.matrix(D),
      identical(dim(D), dim(Y)),
      all(is.finite(D))
    )

    U <- cbind(
      Ub, Ug,
      phi = rowSums(D)
    )
  }

  storage.mode(U) <- "double"

  U
}


# ============================================================
# Compound Poisson-Gamma family specification
# ============================================================

cp_log <- function(power = 1.5) {

  stopifnot(
    length(power) == 1L,
    is.finite(power),
    power > 1,
    power < 2
  )

  list(
    family = "Compound Poisson-Gamma",
    link = "log",

    power = power,

    linkinv = function(eta) {
      exp(eta)
    },

    simulate = function(mu, dispersion) {

      phi <- dispersion

      stopifnot(
        all(is.finite(mu)),
        all(mu > 0),
        length(phi) == 1L,
        is.finite(phi),
        phi > 0
      )

      lambda <- mu^(2 - power) /
        (phi * (2 - power))

      shape <- (2 - power) / (power - 1)

      scale <- phi * (power - 1) *
        mu^(power - 1)

      N <- rpois(length(mu), lambda)

      y <- numeric(length(mu))

      positive <- which(N > 0)

      if (length(positive) > 0L) {
        y[positive] <- rgamma(
          length(positive),
          shape = N[positive] * shape,
          scale = scale[positive]
        )
      }

      y
    },

    derivatives = function(y, mu, dispersion) {

      phi <- dispersion

      g <- (y - mu) *
        mu^(1 - power) / phi

      h <- (
        (1 - power) * y * mu^(1 - power) -
        (2 - power) * mu^(2 - power)
      ) / phi

      list(
        g = g,
        h = h,
        u_dispersion = NULL
      )
    },

    score_batch = function(
        Y, X, Z, mu, theta
    ) {
      cp_score_batch(
        Y = Y,
        X = X,
        Z = Z,
        mu = mu,
        phi = theta,
        power = power
      )
    },

    dispersion = TRUE,
    dispersion_name = "phi"
  )
}


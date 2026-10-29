
# ============================================================
# COMPOUND POISSON-GAMMA FAMILY
# Log link; 1 < power < 2
#
# Family-specific parameters:
#   phi   > 0
#   power in (1, 2)
#
# Both nuisance scores are always included.
# All covariance scores are evaluated at Gamma = 0.
# ============================================================


# ============================================================
# 1. Family parameter validation
# ============================================================

cp_validate_params <- function(phi, power) {

  stopifnot(
    is.numeric(phi),
    length(phi) == 1L,
    is.finite(phi),
    phi > 0,

    is.numeric(power),
    length(power) == 1L,
    is.finite(power),
    power > 1,
    power < 2
  )

  invisible(TRUE)
}


# ============================================================
# 2. Compound Poisson-Gamma simulation
# ============================================================

cp_simulate <- function(mu, phi, power) {

  cp_validate_params(phi, power)

  stopifnot(
    is.numeric(mu),
    all(is.finite(mu)),
    all(mu > 0)
  )

  lambda <- mu^(2 - power) /
    (phi * (2 - power))

  shape <- (2 - power) / (power - 1)

  scale <- phi * (power - 1) *
    mu^(power - 1)

  if (any(!is.finite(lambda)) ||
      any(!is.finite(scale))) {
    stop("Non-finite compound Poisson parameters.")
  }

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
}


# ============================================================
# 3. Dispersion and power scores
#
# Uses the latent Poisson-count representation:
#
#   U_psi(y) = E[U_psi(y, N) | Y = y]
#
# For positive y, conditional expectations are
# calculated using an adaptively truncated series.
# ============================================================

cp_nuisance_score <- function(
    y, mu, phi, power,
    tol = 1e-12,
    max_terms = 100000L
) {

  cp_validate_params(phi, power)

  stopifnot(
    is.numeric(y),
    all(is.finite(y)),
    all(y >= 0),
    is.numeric(mu),
    all(is.finite(mu)),
    all(mu > 0),
    is.numeric(tol),
    length(tol) == 1L,
    is.finite(tol),
    tol > 0,
    tol < 1,
    length(max_terms) == 1L,
    is.finite(max_terms),
    max_terms >= 30,
    max_terms == as.integer(max_terms)
  )

  len <- max(length(y), length(mu))

  stopifnot(
    length(y) %in% c(1L, len),
    length(mu) %in% c(1L, len)
  )

  y <- rep_len(y, len)
  mu <- rep_len(mu, len)

  alpha <- (2 - power) / (power - 1)

  alpha_prime <- -1 / (power - 1)^2

  lambda <- mu^(2 - power) /
    (phi * (2 - power))

  scale <- phi * (power - 1) *
    mu^(power - 1)

  Lp <- -log(mu) + 1 / (2 - power)

  Sp <- log(mu) + 1 / (power - 1)

  Uphi <- numeric(len)
  Up <- numeric(len)

  # ----------------------------------------------------------
  # Exact scores at y = 0
  # ----------------------------------------------------------

  zero <- y == 0

  Uphi[zero] <- lambda[zero] / phi
  Up[zero] <- -lambda[zero] * Lp[zero]

  # ----------------------------------------------------------
  # Positive observations
  # ----------------------------------------------------------

  positive <- which(!zero)

  for (i in positive) {

    yi <- y[i]
    li <- lambda[i]
    si <- scale[i]

    log_y_s <- log(yi) - log(si)

    log_weight <- function(n) {
      n * log(li) -
        lgamma(n + 1) -
        lgamma(n * alpha) +
        (n * alpha - 1) * log(yi) -
        n * alpha * log(si)
    }

    logw <- numeric(0)
    peak <- -Inf
    below <- 0L
    converged <- FALSE

    min_terms <- max(30, ceiling(li))

    if (!is.finite(min_terms) ||
        min_terms > max_terms) {
      stop("Series range exceeds max_terms.")
    }

    for (n in seq_len(max_terms)) {

      lw <- log_weight(n)

      if (!is.finite(lw)) {
        stop("Non-finite series log weight.")
      }

      logw[n] <- lw

      if (lw > peak) {
        peak <- lw
        below <- 0L
      } else if (lw < peak + log(tol)) {
        below <- below + 1L
      } else {
        below <- 0L
      }

      if (n >= min_terms && below >= 30L) {
        converged <- TRUE
        break
      }
    }

    if (!converged) {
      stop(
        "CP series failed to converge; ",
        "increase max_terms or inspect parameters."
      )
    }

    n <- seq_along(logw)

    # Normalize conditional weights
    shift <- max(logw)

    w <- exp(logw - shift)
    w <- w / sum(w)

    # Dispersion score conditional on N
    score_phi_n <- (
      li - n - n * alpha + yi / si
    ) / phi

    # Power score conditional on N
    score_p_n <-
      (n - li) * Lp[i] +
      n * alpha_prime *
        (log_y_s - digamma(n * alpha)) +
      (yi / si - n * alpha) * Sp[i]

    Uphi[i] <- sum(w * score_phi_n)
    Up[i] <- sum(w * score_p_n)
  }

  result <- cbind(
    phi = Uphi,
    power = Up
  )

  if (any(!is.finite(result))) {
    stop("Non-finite compound Poisson-Gamma scores.")
  }

  result
}


# ============================================================
# 4. Observation-level derivatives
# ============================================================

cp_derivatives <- function(
    y, mu, phi, power
) {

  cp_validate_params(phi, power)

  stopifnot(
    is.numeric(y),
    is.numeric(mu),
    length(y) == length(mu),
    all(is.finite(y)),
    all(y >= 0),
    all(is.finite(mu)),
    all(mu > 0)
  )

  # First derivative with respect to eta
  g <- (y - mu) *
    mu^(1 - power) / phi

  # Second derivative with respect to eta
  h <- (
    (1 - power) * y * mu^(1 - power) -
    (2 - power) * mu^(2 - power)
  ) / phi

  nuisance <- cp_nuisance_score(
    y = y,
    mu = mu,
    phi = phi,
    power = power
  )

  list(
    g = g,
    h = h,
    u_family = nuisance
  )
}


# ============================================================
# 5. Vectorized cluster score at Gamma = 0
# ============================================================

cp_score_batch <- function(
    Y, X, Z, mu, phi, power
) {

  cp_validate_params(phi, power)

  stopifnot(
    is.matrix(Y),
    is.numeric(Y),
    all(is.finite(Y)),
    all(Y >= 0),

    is.matrix(X),
    is.numeric(X),
    all(is.finite(X)),

    is.matrix(Z),
    is.numeric(Z),
    all(is.finite(Z)),

    ncol(Y) == nrow(X),
    nrow(X) == nrow(Z),

    length(mu) == ncol(Y),
    all(is.finite(mu)),
    all(mu > 0)
  )

  B <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(X)
  q <- ncol(Z)

  # ----------------------------------------------------------
  # First derivative with respect to eta
  # ----------------------------------------------------------

  G <- sweep(Y, 2L, mu, "-")

  G <- sweep(
    G, 2L,
    mu^(1 - power) / phi,
    "*"
  )

  # ----------------------------------------------------------
  # Second derivative with respect to eta
  # ----------------------------------------------------------

  H <- sweep(
    Y, 2L,
    (1 - power) * mu^(1 - power),
    "*"
  )

  H <- sweep(
    H, 2L,
    (2 - power) * mu^(2 - power),
    "-"
  )

  H <- H / phi

  # ----------------------------------------------------------
  # Fixed-effect scores
  # ----------------------------------------------------------

  Ub <- G %*% X

  # ----------------------------------------------------------
  # Covariance parameter scores
  # ----------------------------------------------------------

  S <- G %*% Z

  n_gamma <- q * (q + 1L) / 2L

  Ug <- matrix(
    0,
    nrow = B,
    ncol = n_gamma
  )

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

  # ----------------------------------------------------------
  # Family-specific scores: phi and power
  # ----------------------------------------------------------

  nuisance <- cp_nuisance_score(
    y = as.vector(Y),
    mu = rep(mu, each = B),
    phi = phi,
    power = power
  )

  Uphi <- rowSums(
    matrix(
      nuisance[, "phi"],
      nrow = B,
      ncol = m
    )
  )

  Up <- rowSums(
    matrix(
      nuisance[, "power"],
      nrow = B,
      ncol = m
    )
  )

  # ----------------------------------------------------------
  # Assemble complete score matrix
  # ----------------------------------------------------------

  beta_names <- colnames(X)

  if (is.null(beta_names)) {
    beta_names <- paste0(
      "beta", seq_len(p)
    )
  }

  colnames(Ub) <- beta_names
  colnames(Ug) <- gamma_names

  U <- cbind(
    Ub,
    Ug,
    phi = Uphi,
    power = Up
  )

  storage.mode(U) <- "double"

  U
}


# ============================================================
# 6. Compound Poisson-Gamma family specification
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

    family_param_names = c("phi", "power"),

    linkinv = function(eta) {
      exp(eta)
    },

    # Simulate one cluster at Gamma = 0
    simulate = function(mu, dispersion) {
      cp_simulate(
        mu = mu,
        phi = dispersion,
        power = power
      )
    },

    # Observation-level derivatives
    derivatives = function(
        y, mu, dispersion
    ) {
      cp_derivatives(
        y = y,
        mu = mu,
        phi = dispersion,
        power = power
      )
    },

    # Vectorized Monte Carlo score
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

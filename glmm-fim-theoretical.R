# Deterministic Fisher information for a negative-binomial log-link GLMM
# at the random-effects covariance boundary Gamma = 0.
# X, Z: lists of cluster-specific model matrices.

.nb_theta_info <- function(mu, theta, tail_tol, max_terms) {
  # I(theta,theta) = sum_{k>=0} P(Y>k)/(theta+k)^2 - 1/theta + 1/(theta+mu)
  K <- min(64L, max_terms)
  repeat {
    k <- seq.int(0L, K - 1L)
    surv <- pnbinom(k, size = theta, mu = mu, lower.tail = FALSE)
    value <- sum(surv / (theta + k)^2) - 1/theta + 1/(theta + mu)
    # Bound omitted positive series by P(Y > K-1) * trigamma(theta+K).
    bound <- surv[K] * trigamma(theta + K)
    if (bound <= tail_tol) return(value)
    if (K >= max_terms) stop("theta information series failed to converge")
    K <- min(max_terms, 2L * K)
  }
}

#' Theoretical boundary Fisher information for an NB-log GLMM
#'
#' @param X,Z Equal-length lists of cluster-specific model matrices.
#' @param beta Fixed-effect coefficients, in the column order of X.
#' @param family Currently only "nbinom" or "nbinom_log".
#' @param dispersion Positive negative-binomial size parameter theta.
#' @param tail_tol Absolute upper bound on omitted series terms per observation.
#' @param max_terms Maximum number of terms per distinct mean.
#' @return List containing I_full, I_metric, and parameter indices.
#' @export
glmm_fim_theoretical <- function(X, Z, beta, family = "nbinom",
                                 dispersion, tail_tol = 1e-12,
                                 max_terms = 1000000L) {
  if (!is.character(family) || length(family) != 1L ||
      !family %in% c("nbinom", "nbinom_log")) {
    stop("Only family = 'nbinom' (log link) is currently supported")
  }
  theta <- dispersion
  if (!is.numeric(theta) || length(theta) != 1L ||
      !is.finite(theta) || theta <= 0) stop("dispersion must be positive")
  if (!is.list(X) || !is.list(Z) || !length(X) || length(X) != length(Z))
    stop("X and Z must be nonempty lists of equal length")
  if (!is.matrix(X[[1L]]) || !is.matrix(Z[[1L]]))
    stop("X and Z must contain matrices")
  p <- ncol(X[[1L]]); q <- ncol(Z[[1L]])
  if (!p || !q || length(beta) != p || any(!is.finite(beta)))
    stop("Invalid dimensions or beta")
  if (!is.numeric(tail_tol) || length(tail_tol) != 1L ||
      !is.finite(tail_tol) || tail_tol <= 0 || tail_tol >= 1)
    stop("tail_tol must be between 0 and 1")
  if (length(max_terms) != 1L || !is.finite(max_terms) ||
      max_terms < 1 || max_terms != as.integer(max_terms))
    stop("max_terms must be a positive integer")
  max_terms <- as.integer(max_terms)

  pairs <- which(upper.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
  ng <- nrow(pairs)
  factors <- ifelse(pairs[, 1L] == pairs[, 2L], 0.5, 1)
  bn <- colnames(X[[1L]])
  if (is.null(bn)) bn <- paste0("b", seq_len(p) - 1L)
  gn <- paste0("g", pairs[, 1L], pairs[, 2L])
  nm <- c(bn, gn, "theta")
  bi <- seq_len(p); gi <- p + seq_len(ng); ti <- p + ng + 1L
  I <- matrix(0, ti, ti, dimnames = list(nm, nm))
  cache <- new.env(parent = emptyenv())

  for (i in seq_along(X)) {
    Xi <- X[[i]]; Zi <- Z[[i]]
    if (!is.matrix(Xi) || !is.matrix(Zi) || ncol(Xi) != p ||
        ncol(Zi) != q || nrow(Xi) != nrow(Zi) || !nrow(Xi) ||
        any(!is.finite(Xi)) || any(!is.finite(Zi)))
      stop("Invalid matrices in cluster ", i)
    mu <- as.vector(exp(Xi %*% beta))
    if (any(!is.finite(mu)) || any(mu <= 0)) stop("Non-finite mean")
    v <- theta * mu / (theta + mu)
    kappa <- theta * mu * (2 * theta * mu + 3 * mu + theta) /
      (theta + mu)^2
    w <- mu^2 / (theta + mu)^2
    M <- crossprod(Zi, v * Zi)
    I[bi, bi] <- I[bi, bi] + crossprod(Xi, v * Xi)
    for (u in seq_len(ng)) {
      a <- pairs[u, 1L]; b <- pairs[u, 2L]
      zab <- Zi[, a] * Zi[, b]
      I[bi, gi[u]] <- I[bi, gi[u]] +
        factors[u] * as.vector(crossprod(Xi, v * zab))
      I[gi[u], ti] <- I[gi[u], ti] - factors[u] * sum(w * zab)
      for (h in seq_len(ng)) {
        c <- pairs[h, 1L]; d <- pairs[h, 2L]
        I[gi[u], gi[h]] <- I[gi[u], gi[h]] + factors[u] * factors[h] * (
          M[a, c] * M[b, d] + M[a, d] * M[b, c] +
          sum((kappa - 2 * v^2) * zab * Zi[, c] * Zi[, d])
        )
      }
    }
    for (m in mu) {
      key <- sprintf("%.17g", m)
      if (!exists(key, envir = cache, inherits = FALSE))
        assign(key, .nb_theta_info(m, theta, tail_tol, max_terms), envir = cache)
      I[ti, ti] <- I[ti, ti] + get(key, envir = cache, inherits = FALSE)
    }
  }
  I[gi, bi] <- t(I[bi, gi, drop = FALSE])
  I[ti, gi] <- I[gi, ti]
  # Schur complement removes beta and theta as nuisance parameters.
  nuisance <- c(bi, ti)
  Inuis <- I[nuisance, nuisance, drop = FALSE]
  Ig_nuis <- I[gi, nuisance, drop = FALSE]
  metric <- I[gi, gi, drop = FALSE] -
    Ig_nuis %*% solve(Inuis, t(Ig_nuis))
  metric <- (metric + t(metric)) / 2
  list(I_full = I, I_metric = metric,
       beta_idx = bi, gamma_idx = gi, dispersion_idx = ti,
       method = "theoretical", boundary = TRUE)
}

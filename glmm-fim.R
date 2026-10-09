
glmm_fim <- function(
    X,
    Z,
    beta,
    family,
    dispersion = NULL,
    nsim = 10000,
    seed = NULL,
    engine = c("auto", "batch", "individual"),
    chunk_size = 1000L
) {

  # -----------------------------
  # Basic checks
  # -----------------------------

  engine <- match.arg(engine)

  if (!is.list(X) || !is.list(Z)) {
    stop("X and Z must be lists of cluster-specific design matrices.")
  }

  n_units <- length(X)

  if (length(Z) != n_units) {
    stop("X and Z must contain the same number of clusters.")
  }

  if (n_units == 0L) {
    stop("At least one cluster is required.")
  }

  if (length(beta) != ncol(X[[1]])) {
    stop("Length of beta must equal the number of columns of X.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # -----------------------------
  # Storage
  # -----------------------------

  I_cluster <- vector("list", n_units)

  score_mean_cluster <- vector("list", n_units)
  score_mcse_cluster <- vector("list", n_units)

  # -----------------------------
  # Monte Carlo by cluster
  # -----------------------------

  for (i in seq_len(n_units)) {

    X_i <- X[[i]]
    Z_i <- Z[[i]]

    # Cluster-specific dimension checks
    if (nrow(X_i) != nrow(Z_i)) {
      stop(
        "X and Z must have the same number of rows within cluster ",
        i,
        "."
      )
    }

    if (ncol(X_i) != length(beta)) {
      stop(
        "All X matrices must have ",
        length(beta),
        " columns."
      )
    }

    if (i > 1L && ncol(Z_i) != ncol(Z[[1]])) {
      stop(
        "All Z matrices must have the same number of columns."
      )
    }

    # Monte Carlo information for cluster i
    res_i <- mc_fisher_cluster(
      X          = X_i,
      Z          = Z_i,
      beta       = beta,
      family     = family,
      dispersion = dispersion,
      nsim       = nsim,
      engine     = engine,
      chunk_size = chunk_size
    )

    I_cluster[[i]] <- res_i$information

    score_mean_cluster[[i]] <- res_i$score_mean
    score_mcse_cluster[[i]] <- res_i$score_mcse
  }

  # -----------------------------
  # Full Fisher information
  # -----------------------------

  I_full <- Reduce("+", I_cluster)

  # -----------------------------
  # Diagnostics for full score
  # -----------------------------

  score_mean_full <-
    Reduce("+", score_mean_cluster)

  # Since clusters are simulated independently,
  # variances of their MC mean estimates add.
  score_mcse_full <- sqrt(
    Reduce(
      "+",
      lapply(
        score_mcse_cluster,
        function(x) x^2
      )
    )
  )

  # Standardized mean-score diagnostic
  standardized_mean_score <-
    score_mean_full / score_mcse_full

  # -----------------------------
  # Parameter information
  # -----------------------------

  parameter_names <- colnames(I_full)

  p <- ncol(X[[1]])
  q <- ncol(Z[[1]])

  n_gamma <- q * (q + 1) / 2

  n_dispersion <-
    if (isTRUE(family$dispersion)) 1L else 0L

  beta_idx <- seq_len(p)

  gamma_idx <- p + seq_len(n_gamma)

  if (n_dispersion > 0L) {
    dispersion_idx <- p + n_gamma + seq_len(n_dispersion)
  } else {
    dispersion_idx <- integer(0)
  }

  # -----------------------------
  # Return glmm_fim object
  # -----------------------------

  structure(
    list(
      I_full          = I_full,
      I_cluster       = I_cluster,

      score_mean      = score_mean_full,
      score_mcse      = score_mcse_full,
      standardized_mean_score = standardized_mean_score,

      score_mean_cluster = score_mean_cluster,
      score_mcse_cluster = score_mcse_cluster,

      parameters      = parameter_names,

      # Parameter blocks
      beta_idx        = beta_idx,
      gamma_idx       = gamma_idx,
      dispersion_idx  = dispersion_idx,

      # Monte Carlo settings
      nsim            = nsim,
      n_units         = n_units,
      engine          = engine,
      chunk_size      = chunk_size,

      family          = family$family,
      link            = family$link,

      p               = p,
      q               = q,
      n_gamma         = n_gamma,
      n_dispersion    = n_dispersion
    ),
    class = "glmm_fim"
  )
}


mc_fisher_cluster <- function(
    X,
    Z,
    beta,
    family,
    dispersion,
    nsim,
    engine = c("auto", "batch", "individual"),
    chunk_size = 1000L
) {

  # --------------------------------------------------
  # 1. Validate arguments
  # --------------------------------------------------

  engine <- match.arg(engine)

  if (length(nsim) != 1L ||
      !is.numeric(nsim) ||
      is.na(nsim) ||
      !is.finite(nsim) ||
      nsim < 2 ||
      nsim != floor(nsim)) {
    stop("nsim must be an integer >= 2.")
  }

  if (length(chunk_size) != 1L ||
      !is.numeric(chunk_size) ||
      is.na(chunk_size) ||
      !is.finite(chunk_size) ||
      chunk_size < 1 ||
      chunk_size != floor(chunk_size)) {
    stop("chunk_size must be a positive integer.")
  }

  if (nsim > .Machine$integer.max ||
      chunk_size > .Machine$integer.max) {
    stop("nsim and chunk_size must fit within integer limits.")
  }

  nsim <- as.integer(nsim)
  chunk_size <- min(as.integer(chunk_size), nsim)

  # --------------------------------------------------
  # 2. Select Monte Carlo engine
  # --------------------------------------------------

  has_batch <- is.function(family$score_batch)

  if (engine == "auto") {
    engine <- if (has_batch) "batch" else "individual"
  }

  if (engine == "batch" && !has_batch) {
    stop(
      "This family does not provide a score_batch() method. ",
      "Use engine = 'individual' or 'auto'."
    )
  }

  # --------------------------------------------------
  # 3. Mean at Gamma = 0
  # --------------------------------------------------

  eta <- drop(X %*% beta)
  mu  <- family$linkinv(eta)

  m <- length(mu)

  # --------------------------------------------------
  # 4. Compute scores for one chunk
  # --------------------------------------------------

  score_chunk <- function(B) {

    if (engine == "batch") {

      # Generate B independent response vectors.
      # Rows = simulations; columns = observations.

      Y <- matrix(
        family$simulate(
          mu = rep(mu, each = B),
          dispersion = dispersion
        ),
        nrow = B,
        ncol = m
      )

      U <- family$score_batch(
        Y = Y,
        X = X,
        Z = Z,
        mu = mu,
        theta = dispersion
      )

    } else {

      # General individual-score implementation

      scores <- vector("list", B)

      for (r in seq_len(B)) {

        y <- family$simulate(
          mu = mu,
          dispersion = dispersion
        )

        scores[[r]] <- score_glmm(
          y = y,
          X = X,
          Z = Z,
          beta = beta,
          family = family,
          dispersion = dispersion
        )
      }

      U <- do.call(rbind, scores)
    }

    # Validate returned score matrix

    if (!is.matrix(U) ||
        nrow(U) != B ||
        !is.numeric(U) ||
        any(!is.finite(U))) {
      stop("Score computation returned an invalid matrix.")
    }

    U
  }

  # --------------------------------------------------
  # 5. Initialize running moments
  # --------------------------------------------------

  N <- 0L
  mean_U <- NULL
  M <- NULL
  parameter_names <- NULL

  # --------------------------------------------------
  # 6. Monte Carlo simulation in chunks
  # --------------------------------------------------

  while (N < nsim) {

    B <- min(chunk_size, nsim - N)

    U <- score_chunk(B)

    if (N == 0L) {

      parameter_names <- colnames(U)

      if (is.null(parameter_names)) {
        stop("Score matrix must have parameter names.")
      }

    } else {

      if (ncol(U) != length(parameter_names) ||
          !identical(colnames(U), parameter_names)) {
        stop("Inconsistent score columns across chunks.")
      }
    }

    # Chunk mean

    mean_B <- colMeans(U)

    # Center scores within the chunk

    U_centered <- sweep(
      U,
      MARGIN = 2L,
      STATS = mean_B,
      FUN = "-"
    )

    # Chunk centered cross-product

    M_B <- crossprod(U_centered)

    # Combine chunk statistics

    if (N == 0L) {

      mean_U <- mean_B
      M <- M_B

    } else {

      delta <- mean_B - mean_U
      total <- N + B

      M <- M + M_B +
        (N * B / total) * tcrossprod(delta)

      mean_U <- mean_U +
        (B / total) * delta
    }

    N <- N + B
  }

  # --------------------------------------------------
  # 7. Fisher information and diagnostics
  # --------------------------------------------------

  information <- M / nsim

  score_mean <- mean_U
  names(score_mean) <- parameter_names

  # Sample SD divided by sqrt(nsim)

  score_mcse <- sqrt(
    diag(M) / (nsim - 1)
  ) / sqrt(nsim)

  names(score_mcse) <- parameter_names

  # --------------------------------------------------
  # 8. Return
  # --------------------------------------------------

  list(
    information = information,
    score_mean  = score_mean,
    score_mcse  = score_mcse
  )
}

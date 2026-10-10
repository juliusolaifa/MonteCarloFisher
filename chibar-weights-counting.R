#' Chi-bar-square weights by Gaussian cone counting (2x2 PSD cone)
#'
#' @param information Efficient 3x3 information matrix for
#'   gamma = (gamma11, gamma12, gamma22), in that order.
#' @param nsim Number of Gaussian draws.
#' @param seed Optional RNG seed; restores the caller's RNG state.
#' @param chunk_size Draws per batch to limit memory use.
#' @return List with w = (w0,w1,w2,w3), w0 and w3 Monte Carlo SEs.
#' @export
chibar_weights_counting <- function(information, nsim = 100000L,
                                   seed = NULL, chunk_size = 100000L) {
  I <- information
  if (!is.matrix(I) || !identical(dim(I), c(3L, 3L)) ||
      any(!is.finite(I)) || max(abs(I - t(I))) > 1e-8)
    stop("information must be a finite symmetric 3x3 matrix")
  if (inherits(try(chol(I), silent = TRUE), "try-error"))
    stop("information must be positive definite")
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1 || nsim != floor(nsim))
    stop("nsim must be a positive integer")
  if (length(chunk_size) != 1L || !is.finite(chunk_size) ||
      chunk_size < 1 || chunk_size != floor(chunk_size))
    stop("chunk_size must be a positive integer")
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed)) stop("Invalid seed")
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  # Rows of Z have covariance I^{-1}; rows of U = Z I have covariance I.
  root <- chol(solve(I))
  n0 <- 0; n3 <- 0; remaining <- nsim
  while (remaining > 0) {
    B <- min(remaining, chunk_size)
    Z <- matrix(rnorm(B * 3L), ncol = 3L) %*% root
    inC <- (Z[, 1L] + Z[, 3L] >= 0) &
      (Z[, 1L] * Z[, 3L] - Z[, 2L]^2 >= 0)
    U <- Z %*% I
    inP <- (U[, 1L] + U[, 3L] <= 0) &
      (U[, 1L] * U[, 3L] - (U[, 2L]/2)^2 >= 0)
    n3 <- n3 + sum(inC)
    n0 <- n0 + sum(inP)
    remaining <- remaining - B
  }
  w0 <- n0/nsim; w3 <- n3/nsim
  w <- c(w0 = w0, w1 = 0.5 - w3, w2 = 0.5 - w0, w3 = w3)
  list(w = w, se = c(w0 = sqrt(w0 * (1-w0)/nsim),
                      w3 = sqrt(w3 * (1-w3)/nsim)),
       nsim = nsim, method = "counting")
}

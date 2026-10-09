glmm_fim <- function(
    fixed,
    random,
    data,
    cluster,
    beta,
    family,
    dispersion = NULL,
    nsim = 10000,
    seed = NULL,
    engine = c("auto", "batch", "individual"),
    chunk_size = 1000L
) {

  # ---------------------------------
  # 1. Validate inputs
  # ---------------------------------

  engine <- match.arg(engine)

  if (!inherits(fixed, "formula")) {
    stop("'fixed' must be a formula.")
  }

  if (!inherits(random, "formula")) {
    stop("'random' must be a formula.")
  }

  if (length(fixed) != 2L || length(random) != 2L) {
    stop("'fixed' and 'random' must be one-sided formulas.")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  if (!is.character(cluster) ||
      length(cluster) != 1L ||
      is.na(cluster) ||
      !cluster %in% names(data)) {
    stop("'cluster' must name a column in 'data'.")
  }

  if (nrow(data) == 0L) {
    stop("'data' must contain at least one observation.")
  }

  if (anyNA(data[[cluster]])) {
    stop("Cluster identifiers cannot contain missing values.")
  }

  # ---------------------------------
  # 2. Construct design matrices
  # ---------------------------------

  # model.matrix() handles intercepts,
  # factors, interactions, and contrasts.

  X <- model.matrix(
    object = fixed,
    data   = data,
    na.action = na.fail
  )

  Z <- model.matrix(
    object = random,
    data   = data,
    na.action = na.fail
  )

  # ---------------------------------
  # 3. Validate matrix dimensions
  # ---------------------------------

  if (length(beta) != ncol(X)) {
    stop(
      "Length of 'beta' must equal the ",
      "number of columns in the fixed-effects ",
      "design matrix (", ncol(X), ")."
    )
  }

  if (ncol(Z) == 0L) {
    stop("The random-effects formula must produce at least one column.")
  }

  # ---------------------------------
  # 4. Split matrices by cluster
  # ---------------------------------

  cluster_id <- data[[cluster]]

  idx <- split(
    seq_len(nrow(data)),
    factor(cluster_id, levels = unique(cluster_id)),
    drop = TRUE
  )

  X_list <- lapply(
    idx,
    function(i) X[i, , drop = FALSE]
  )

  Z_list <- lapply(
    idx,
    function(i) Z[i, , drop = FALSE]
  )

  # ---------------------------------
  # 5. Resolve distribution family
  # ---------------------------------

  if (is.character(family)) {

    if (length(family) != 1L || is.na(family)) {
      stop("'family' must be a single family name.")
    }

    family <- switch(
      family,
      nbinom = nbinom_log(),
      stop("Unsupported family: ", family)
    )
  }

  # ---------------------------------
  # 6. Call existing matrix engine
  # ---------------------------------

  glmm_fim_matrix(
    X          = X_list,
    Z          = Z_list,
    beta       = beta,
    family     = family,
    dispersion = dispersion,
    nsim       = nsim,
    seed       = seed,
    engine     = engine,
    chunk_size = chunk_size
  )
}

check_fim <- function(
    object,
    score_tol = 3,
    eigen_tol = 1e-8,
    digits = 3
) {
  
  if (!inherits(object, "glmm_fim")) {
    stop("object must inherit from class 'glmm_fim'.")
  }
  
  # --------------------------------------------------
  # 1. Full-score diagnostic
  # --------------------------------------------------
  
  z_full <- object$standardized_mean_score
  
  max_z_full <- max(
    abs(z_full),
    na.rm = TRUE
  )
  
  full_score_ok <- max_z_full <= score_tol
  
  
  # --------------------------------------------------
  # 2. Cluster-level score diagnostics
  # --------------------------------------------------
  
  z_cluster <- Map(
    function(mean_i, mcse_i) {
      mean_i / mcse_i
    },
    object$score_mean_cluster,
    object$score_mcse_cluster
  )
  
  max_z_cluster <- vapply(
    z_cluster,
    function(z) {
      max(abs(z), na.rm = TRUE)
    },
    numeric(1)
  )
  
  worst_cluster <- which.max(max_z_cluster)
  
  worst_cluster_z <- max_z_cluster[worst_cluster]
  
  cluster_score_ok <- worst_cluster_z <= score_tol
  
  
  # --------------------------------------------------
  # 3. Symmetry diagnostic
  # --------------------------------------------------
  
  symmetry_error <- max(
    abs(object$I_full - t(object$I_full))
  )
  
  symmetry_ok <- symmetry_error < 1e-10
  
  
  # --------------------------------------------------
  # 4. Positive-semidefinite diagnostic
  # --------------------------------------------------
  
  eig <- eigen(
    object$I_full,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  min_eigenvalue <- min(eig)
  
  psd_ok <- min_eigenvalue >= -eigen_tol
  
  
  # --------------------------------------------------
  # Print
  # --------------------------------------------------
  
  cat("GLMM FIM Monte Carlo Diagnostics\n")
  cat("--------------------------------\n")
  
  cat(
    "Full-score maximum |z|:       ",
    round(max_z_full, 2),
    "  ",
    if (full_score_ok) "[OK]" else "[CHECK]",
    "\n",
    sep = ""
  )
  
  cat(
    "Worst cluster maximum |z|:    ",
    round(worst_cluster_z, 2),
    "  ",
    if (cluster_score_ok) "[OK]" else "[CHECK]",
    "\n",
    sep = ""
  )
  
  cat(
    "Worst cluster:                ",
    worst_cluster,
    "\n",
    sep = ""
  )
  
  cat(
    "Maximum symmetry error:       ",
    signif(symmetry_error, digits),
    "  ",
    if (symmetry_ok) "[OK]" else "[CHECK]",
    "\n",
    sep = ""
  )
  
  cat(
    "Minimum eigenvalue:           ",
    signif(min_eigenvalue, digits),
    "  ",
    if (psd_ok) "[OK]" else "[CHECK]",
    "\n",
    sep = ""
  )
  
  
  # --------------------------------------------------
  # Return diagnostics invisibly
  # --------------------------------------------------
  
  invisible(
    list(
      full = list(
        standardized_mean_score = z_full,
        max_abs_z = max_z_full,
        ok = full_score_ok
      ),
      
      cluster = list(
        standardized_mean_score = z_cluster,
        max_abs_z = max_z_cluster,
        worst_cluster = worst_cluster,
        worst_cluster_z = worst_cluster_z,
        ok = cluster_score_ok
      ),
      
      symmetry = list(
        error = symmetry_error,
        ok = symmetry_ok
      ),
      
      eigen = list(
        values = eig,
        minimum = min_eigenvalue,
        ok = psd_ok
      )
    )
  )
}

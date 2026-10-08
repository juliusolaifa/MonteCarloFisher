print.glmm_fim <- function(x, digits = 3, ...) {
  
  cat("GLMM Fisher Information\n")
  cat("-----------------------\n")
  
  cat("Family:                ", x$family, "\n", sep = "")
  cat("Link:                  ", x$link, "\n", sep = "")
  cat("Fixed effects:         ", x$p, "\n", sep = "")
  cat("Random effects:        ", x$q, "\n", sep = "")
  cat("Covariance parameters: ", x$n_gamma, "\n", sep = "")
  cat("Dispersion parameters: ", x$n_dispersion, "\n", sep = "")
  cat("Monte Carlo replicates: ", x$nsim, "\n", sep = "")
  cat("Independent units:     ", x$n_units, "\n", sep = "")
  
  cat("\nEstimated Fisher information:\n\n")
  
  print(
    round(x$I_full, digits)
  )
  
  max_score <- max(
    abs(x$standardized_mean_score),
    na.rm = TRUE
  )
  
  cat(
    "\nMaximum standardized mean score: ",
    round(max_score, 2),
    "\n",
    sep = ""
  )
  
  invisible(x)
}

summary.glmm_fim <- function(object, digits = 3, ...) {
  
  diagnostics <- data.frame(
    parameter = object$parameters,
    mean_score = as.numeric(object$score_mean),
    mcse = as.numeric(object$score_mcse),
    standardized = as.numeric(
      object$standardized_mean_score
    ),
    row.names = NULL
  )
  
  cat("Summary: GLMM Fisher Information\n")
  cat("--------------------------------\n")
  
  cat("Family:                ", object$family, "\n", sep = "")
  cat("Link:                  ", object$link, "\n", sep = "")
  cat("Fixed effects:         ", object$p, "\n", sep = "")
  cat("Random effects:        ", object$q, "\n", sep = "")
  cat("Covariance parameters: ", object$n_gamma, "\n", sep = "")
  cat("Dispersion parameters: ", object$n_dispersion, "\n", sep = "")
  cat("Monte Carlo replicates: ", object$nsim, "\n", sep = "")
  cat("Independent units:     ", object$n_units, "\n", sep = "")
  
  cat("\nEstimated Fisher information:\n\n")
  
  print(
    round(object$I_full, digits)
  )
  
  cat("\nMonte Carlo score diagnostics:\n\n")
  
  diagnostics_print <- diagnostics
  
  diagnostics_print$mean_score <-
    signif(diagnostics_print$mean_score, digits)
  
  diagnostics_print$mcse <-
    signif(diagnostics_print$mcse, digits)
  
  diagnostics_print$standardized <-
    round(diagnostics_print$standardized, 2)
  
  print(
    diagnostics_print,
    row.names = FALSE
  )
  
  invisible(object)
}

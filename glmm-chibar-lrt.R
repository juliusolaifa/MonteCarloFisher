glmm_chibar_lrt <- function(
    logLik_null,
    logLik_alt,
    weights,
    tol = 1e-8
) {
    # Validate log-likelihoods
    if (
        length(logLik_null) != 1L ||
        length(logLik_alt) != 1L ||
        !is.finite(logLik_null) ||
        !is.finite(logLik_alt)
    ) {
        stop("Log-likelihoods must be finite scalars.")
    }

    # Validate weights
    if (
        !is.numeric(weights) ||
        length(weights) < 2L ||
        any(!is.finite(weights)) ||
        any(weights < 0) ||
        abs(sum(weights) - 1) > tol
    ) {
        stop("Invalid chi-bar-square weights.")
    }

    # Likelihood ratio statistic
    statistic <- 2 * (
        as.numeric(logLik_alt) -
        as.numeric(logLik_null)
    )

    # Allow only small numerical discrepancies
    if (statistic < -tol) {
        stop(
            "Alternative log-likelihood is smaller ",
            "than null log-likelihood."
        )
    }

    statistic <- max(0, statistic)

    # Degrees of freedom: 0, 1, ..., K
    df <- seq_along(weights) - 1L

    # Chi-bar-square upper-tail probability
    if (statistic == 0) {
        p_value <- 1
    } else {
        p_value <- sum(
            weights[-1L] *
                pchisq(
                    statistic,
                    df = df[-1L],
                    lower.tail = FALSE
                )
        )
    }

    list(
        statistic = statistic,
        df = df,
        weights = weights,
        p_value = p_value
    )
}

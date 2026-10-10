glmm_chibar_lrt <- function(
    logLik_null,
    logLik_alt,
    weights
) {
    stopifnot(
        length(logLik_null) == 1L,
        length(logLik_alt) == 1L,
        all(is.finite(c(logLik_null, logLik_alt))),
        is.numeric(weights),
        length(weights) >= 2L,
        all(is.finite(weights)),
        all(weights >= 0),
        abs(sum(weights) - 1) < 1e-8
    )

    statistic <- max(
        0,
        2 * (logLik_alt - logLik_null)
    )

    df <- seq_along(weights) - 1L

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
        weights = weights,
        p_value = p_value
    )
}

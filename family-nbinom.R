nbinom_log <- function() {
  
  list(
    family = "Negative binomial",
    link   = "log",
    
    # mu = g^{-1}(eta)
    linkinv = function(eta) {
      exp(eta)
    },
    
    # Simulate one cluster conditional on Gamma = 0
    simulate = function(mu, dispersion) {
      rnbinom(
        n    = length(mu),
        size = dispersion,
        mu   = mu
      )
    },
    
    # Observation-level derivatives
    derivatives = function(y, mu, dispersion) {
      
      theta <- dispersion
      
      # d log f / d eta
      g <- theta * (y - mu) / (theta + mu)
      
      # d^2 log f / d eta^2
      h <- -theta * mu * (theta + y) /
        (theta + mu)^2
      
      # d log f / d theta
      u_dispersion <-
        digamma(y + theta) -
        digamma(theta) +
        log(theta) + 1 -
        log(theta + mu) -
        (y + theta) / (theta + mu)
      
      list(
        g            = g,
        h            = h,
        u_dispersion = u_dispersion
      )
    },
    
    dispersion      = TRUE,
    dispersion_name = "theta"
  )
}

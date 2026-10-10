#' Deterministic chi-bar-square weights for the 2x2 PSD cone
#'
#' Uses 2D spherical quadrature of the angular Gaussian density.
#' This is deterministic numerical integration, not a closed-form formula.
#' @param information Efficient 3x3 information matrix for
#'   gamma = (gamma11, gamma12, gamma22), in that order.
#' @param rel.tol Relative tolerance for numerical integration.
#' @param subdivisions Maximum subdivisions in each 1D integration.
#' @return List with w = (w0,w1,w2,w3) and integration error estimates.
#' @export
chibar_weights_theoretical <- function(information, rel.tol = 1e-7,
                                      subdivisions = 200L) {
  I <- information
  if (!is.matrix(I) || !identical(dim(I), c(3L, 3L)) ||
      any(!is.finite(I)) || max(abs(I - t(I))) > 1e-8)
    stop("information must be a finite symmetric 3x3 matrix")
  if (inherits(try(chol(I), silent = TRUE), "try-error"))
    stop("information must be positive definite")
  if (!is.finite(rel.tol) || length(rel.tol) != 1L ||
      rel.tol <= 0 || rel.tol >= 1) stop("Invalid rel.tol")
  if (length(subdivisions) != 1L || !is.finite(subdivisions) ||
      subdivisions < 1) stop("Invalid subdivisions")

  s <- 1/sqrt(2)
  # Map gamma=(a,b,c) to Lorentz coordinates
  # (t,u,v)=((a+c)/sqrt(2),(a-c)/sqrt(2),sqrt(2)*b).
  A_primal <- rbind(c(s, 0, s), c(s, 0, -s), c(0, sqrt(2), 0))
  # For U = Z I, the polar cone is -dual(PSD); its Lorentz map is
  # (-(U1+U3)/sqrt(2),(U1-U3)/sqrt(2),U2/sqrt(2)).
  A_polar <- rbind(c(-s, 0, -s), c(s, 0, -s), c(0, s, 0))

  cone_probability <- function(Sigma) {
    P <- solve(Sigma)
    normalizer <- 1/(4*pi*sqrt(det(Sigma)))
    # Unit sphere: t=cos(alpha), u=sin(alpha)cos(phi),
    # v=sin(alpha)sin(phi). Lorentz cone: 0 <= alpha <= pi/4.
    inner <- function(alpha) {
      vapply(alpha, function(a) {
        sa <- sin(a); ca <- cos(a)
        f <- function(phi) {
          u <- sa*cos(phi); v <- sa*sin(phi)
          quad <- P[1,1]*ca^2 + P[2,2]*u^2 + P[3,3]*v^2 +
            2*P[1,2]*ca*u + 2*P[1,3]*ca*v + 2*P[2,3]*u*v
          normalizer * quad^(-1.5) * sa
        }
        integrate(f, 0, 2*pi, rel.tol = rel.tol,
                  subdivisions = subdivisions)$value
      }, numeric(1))
    }
    integrate(inner, 0, pi/4, rel.tol = rel.tol,
              subdivisions = subdivisions)
  }
  # Z ~ N(0, I^{-1}); U=ZI ~ N(0,I).
  primal <- cone_probability(A_primal %*% solve(I) %*% t(A_primal))
  polar <- cone_probability(A_polar %*% I %*% t(A_polar))
  w3 <- primal$value; w0 <- polar$value
  w <- c(w0 = w0, w1 = 0.5-w3, w2 = 0.5-w0, w3 = w3)
  if (any(w < -1e-7) || any(w > 1+1e-7))
    warning("Computed weights outside [0,1]; check integration accuracy")
  list(w = w, abs.error = c(w0 = polar$abs.error, w3 = primal$abs.error),
       method = "theoretical_quadrature")
}


#' Get expected openness and 1/openness
#'
#' Get expected openness and 1/openness through the closed-form GIG expectations and `E_openness`.
#'
#' @details
#' The closed-form generalized inverse Gaussian (GIG) expectations (`egig_x` /
#' `egig_inv_x`, evaluated with base R `besselK`) and `E_openness` work for
#' different sets of parameters. To reduce failures in obtaining expected values,
#' first use the closed-form expectations, next perform `E_openness` if they
#' return non-finite values.
#'
#' @param lambda numeric vector of GIG parameter lambda.
#' @param chi numeric vector of GIG parameter chi. Must be positive.
#' @param psi numeric vector of GIG parameter psi. Must be positive.
#' @param alpha numeric vector of inverse-gamma shape parameter.
#'
#' @return A list object of the following expected values:
#'   \describe{
#'     \item{o}{numeric vector of expected opennes.}
#'     \item{inv_o}{numeric vector of expected 1/openness.}
#'   }
#'
#' @export
get_openness <- function(lambda, chi, psi, alpha){
  # Parameters for the closed-form GIG expectations
  mx1 <- matrix(c(lambda, chi, psi), nrow=3, byrow=TRUE)
  # Expected x and 1/x of the GIG distribution
  gig_o <- sapply(as.data.frame(mx1), exec_Egig_x)
  gig_inv_o <- sapply(as.data.frame(mx1), exec_Egig_inv_x)

  # Parameters for E_openness
  mx2 <- matrix(c(gig_o, gig_inv_o, lambda, chi, psi, alpha), nrow=6, byrow=TRUE)
  # Execute E_openness
  sapply(as.data.frame(mx2), exec_E_openness)
}

#' Closed-form expectations of the generalized inverse Gaussian distribution
#'
#' Expected value of `x` and `1/x` for a variable `x ~ GIG(lambda, chi, psi)`,
#' evaluated with base R `besselK`. This replaces `ghyp::Egig` so that the
#' package has no dependency on the `ghyp` package.
#'
#' @details
#' For `x ~ GIG(lambda, chi, psi)` with `w = sqrt(chi * psi)` and the Bessel-K
#' ratio `R = K_{lambda+1}(w) / K_{lambda}(w)`:
#' \deqn{E[x]   = sqrt(chi / psi) * R}
#' \deqn{E[1/x] = sqrt(psi / chi) * R - 2 * lambda / chi}
#' The ratio `R` is computed with `expon.scaled = TRUE` so that the shared
#' `exp(-w)` factor cancels, keeping the ratio numerically stable for large `w`.
#'
#' @param lambda numeric GIG parameter lambda.
#' @param chi numeric GIG parameter chi. Must be positive.
#' @param psi numeric GIG parameter psi. Must be positive.
#'
#' @return The expected value of `x` (`egig_x`) or `1/x` (`egig_inv_x`).
egig_x <- function(lambda, chi, psi){
  w <- sqrt(chi * psi)
  r <- besselK(w, lambda + 1, expon.scaled=TRUE) / besselK(w, lambda, expon.scaled=TRUE)
  sqrt(chi / psi) * r
}

#' @rdname egig_x
egig_inv_x <- function(lambda, chi, psi){
  w <- sqrt(chi * psi)
  r <- besselK(w, lambda + 1, expon.scaled=TRUE) / besselK(w, lambda, expon.scaled=TRUE)
  sqrt(psi / chi) * r - 2 * lambda / chi
}

#' Execute the closed-form GIG expectation of `x`
#'
#' Wrapper for `egig_x` / `egig_inv_x` to execute through `apply` function family.
#' Returns `NA` (rather than a non-finite value) so that `get_openness` falls
#' back to `E_openness` when the closed form is not usable.
#'
#' @param v numeric vector of parameters for GIG in the form `c(lambda, chi, psi)`.
#'
#' @return Expected `x` or `1/x` of the GIG distribution
#'   \itemize{
#'     \item `exec_Egig_x` gives expected `x`
#'     \item `exec_Egig_inv_x` gives expected `1/x`
#'   }
exec_Egig_x <- function(v){
  tryCatch({
    val <- egig_x(lambda=v[1], chi=v[2], psi=v[3])
    if(!is.finite(val)) NA else val
  }, error=function(e){
    NA
  })
}

#' @rdname exec_Egig_x
exec_Egig_inv_x <- function(v){
  tryCatch({
    val <- egig_inv_x(lambda=v[1], chi=v[2], psi=v[3])
    if(!is.finite(val)) NA else val
  }, error=function(e){
    NA
  })
}

#' Execute E_openness
#'
#' Wrapper for `E_openness` to execute through `apply` function family.
#'
#' @param v numeric vector in the form `c(o, inv_o, lambda, chi, psi, alpha)`.
#' If both `o` and `inv_o` are valid, do nothing and return the `o` and `inv_o`.
#' If not, `E_openness` is executed with `lambda`, `chi`, `psi`, and `alpha` as parameters.
#'
#' @return A list object of the following expected values:
#'   \describe{
#'     \item{o}{expected opennes.}
#'     \item{inv_o}{expected 1/openness.}
#'   }
exec_E_openness <- function(v){
  if(is.nan(v[1]) | is.nan(v[2]) | is.infinite(v[1]) | is.infinite(v[2]) | is.na(v[1]) | is.na(v[2])){
    E_openness(lambda=v[3], chi=v[4], psi=v[5], alpha=v[6])
  }else{
    list(o=v[1], inv_o=v[2])
  }
}

#' Calculate expected values of openness and 1/openness
#'
#' Expected values of openness and 1/openness using `stats::integrate` with
#' their density function described by the product of gamma and inverse-gamma distributions.
#'
#' @param lambda numeric vector of GIG parameter lambda.
#' @param chi numeric vector of GIG parameter chi. Must be positive.
#' @param psi numeric vector of GIG parameter psi. Must be positive.
#' @param alpha numeric vector of inverse-gamma shape parameter.
#'
#' @return A list object of the following expected values:
#'   \describe{
#'     \item{o}{numeric vector of expected opennes.}
#'     \item{inv_o}{numeric vector of expected 1/openness.}
#'   }
#'
#' @export
#E_openness <- function(a, b, c, d){
E_openness <- function(lambda, chi, psi, alpha){
  S <- NULL; E_o <- NULL; E_inv_o <- NULL;
  lower <- 10^-100
  a <- lambda + alpha + 1 # gamma shape parameter.
  b <- psi / 2            # gamma rate parameter.
  c <- alpha              # inverse-gamma shape parameter.
  d <- chi / 2            # inverse-gamma rate parameter.

  upper <- max(ceiling(3 * a / b), sqrt(a / b * d / c))
  # Scaling factor. The density is scaled for stable calculation in stats::integrate.
  sf <- 10^8
  tryCatch({
    s0 <- (stats::integrate(density_openness, lower, upper, a=a, b=b, c=c, d=d, S=sf))$value * 10
    # Scaled normalization constant
    S <- stats::integrate(density_openness, lower, upper, a=a, b=b, c=c, d=d, S=(sf / s0))
    # Value for expected openness
    E_o <- stats::integrate(weighted_o, lower, upper, a=a, b=b, c=c, d=d, S=(sf^2 / (S$value * s0)))
    # Value for expected 1/openness
    E_inv_o <- stats::integrate(weighted_inv_o, lower, upper, a=a, b=b, c=c, d=d, S=(sf^2 / (S$value * s0)))

    # Errors
    err_S <- abs(S$abs.error / S$value)
    err_E_o <- abs(E_o$abs.error / E_o$value)
    err_E_inv_o <- abs(E_inv_o$abs.error / E_inv_o$value)
    if(err_S > 0.01){ S <- NULL; message("Large error in S: ", err_S); }
    if(err_E_o > 0.01){ E_o <- NULL; message("Large error in E_o: ", err_E_o); }
    if(err_E_inv_o > 0.01){ E_inv_o <- NULL; message("Large error in E_inv_o: ", err_E_inv_o); }

  }, error=function(e){
    message(e, "in S, E_o or E_inv_o.")
    S <- NULL; E_o <- NULL; E_inv_o <- NULL;
  })

  if(!is.null(S) & !is.null(E_o) & !is.null(E_inv_o)){
    return(list(o=(E_o$value / sf), inv_o=(E_inv_o$value / sf)))
  }else{
    return(list(o=NA, inv_o=NA))
  }

}

#' Density function of openness
#'
#' Density function used for calculating expected values of openness and 1/openness.
#' The function is described by the product of gamma and inverse-gamma distributions.
#'
#' @param x numeric vector of quantiles.
#' @param a gamma shape parameter.
#' @param b gamma rate parameter.
#' @param c inverse-gamma shape parameter.
#' @param d inverse-gamma rate parameter.
#' @param S scaling constant
#'
#' @return
#'   \itemize{
#'     \item `density_openness` gives scaled density.
#'     \item `weighted_o` gives density-weighted openness.
#'     \item `weighted_inv_o` gives density-weighted 1/openness.
#'   }
density_openness <- function(x, a, b, c, d, S=1){
  stats::dgamma(x, shape=a, rate=b) * invgamma::dinvgamma(x, shape=c, rate=d) * S
}

#' Function for density-weighted openness
#'
#' @rdname density_openness
weighted_o <- function(x, a, b, c, d, S=1){
  a / b * density_openness(x, a + 1, b, c, d, S)
}

#' Function for density-weighted 1/openness
#'
#' @rdname density_openness
weighted_inv_o <- function(x, a, b, c, d, S=1){
  c / d * density_openness(x, a, b, c + 1, d, S)
}

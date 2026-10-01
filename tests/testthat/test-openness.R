# Independent reference for the generalized inverse Gaussian (GIG) expectations,
# computed by direct numerical integration of the GIG density. This validates the
# closed-form besselK implementation in openness.R without depending on {ghyp}.
#
# For x ~ GIG(lambda, chi, psi) the (unnormalized) density is
#   f(x) proportional to x^(lambda-1) exp(-(chi/x + psi*x)/2),  x > 0.
gig_ref <- function(lambda, chi, psi, inv = FALSE){
  logk <- function(x) (lambda - 1) * log(x) - (chi / x + psi * x) / 2
  # Shift by the log-density at the mode so integrate() works on an O(1) scale.
  mode <- (lambda - 1 + sqrt((lambda - 1)^2 + chi * psi)) / psi
  if(!is.finite(mode) || mode <= 0) mode <- sqrt(chi / psi)
  lp <- logk(mode)
  dens <- function(x) exp(logk(x) - lp)
  norm <- stats::integrate(dens, 0, Inf, rel.tol = 1e-11)$value
  if(inv){
    stats::integrate(function(x) dens(x) / x, 0, Inf, rel.tol = 1e-11)$value / norm
  }else{
    stats::integrate(function(x) x * dens(x), 0, Inf, rel.tol = 1e-11)$value / norm
  }
}


test_that("get_openness matches the GIG moments (besselK closed form)", {
  alpha <- 10

  # Normal case: both the closed form and E_openness are usable.
  lambda <- -5; chi <- 10; psi <- 2
  ref_o <- gig_ref(lambda, chi, psi)
  ref_inv_o <- gig_ref(lambda, chi, psi, inv = TRUE)
  ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
  open <- get_openness(lambda, chi, psi, alpha)
  # Closed form (used by get_openness) reproduces the numerically-integrated moments.
  expect_equal(unlist(open[1, ]), ref_o, ignore_attr = TRUE)
  expect_equal(unlist(open[2, ]), ref_inv_o, ignore_attr = TRUE)
  # E_openness independently agrees with the same moments.
  expect_equal(ls_E_o$o[1], ref_o, ignore_attr = TRUE, tolerance = 1e-4)
  expect_equal(ls_E_o$inv_o[1], ref_inv_o, ignore_attr = TRUE, tolerance = 1e-4)

  # Case where E_openness cannot be integrated but the closed form still works:
  # get_openness must return the closed-form (GIG) moments.
  lambda <- -8; chi <- 0.01; psi <- 2
  ref_o <- gig_ref(lambda, chi, psi)
  ref_inv_o <- gig_ref(lambda, chi, psi, inv = TRUE)
  ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
  open <- get_openness(lambda, chi, psi, alpha)
  expect_true(is.na(ls_E_o$o[1]))
  expect_true(is.na(ls_E_o$inv_o[1]))
  expect_equal(unlist(open[1, ]), ref_o, ignore_attr = TRUE)
  expect_equal(unlist(open[2, ]), ref_inv_o, ignore_attr = TRUE)
})


test_that("get_openness falls back to E_openness when the closed form overflows", {
  alpha <- 10

  # Large lambda makes besselK overflow (non-finite), so get_openness should
  # return exactly the E_openness result instead.
  for(p in list(c(1000, 2000, 2), c(10000, 20000, 2))){
    lambda <- p[1]; chi <- p[2]; psi <- p[3]
    ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
    open <- get_openness(lambda, chi, psi, alpha)
    expect_false(is.na(ls_E_o$o[1]))
    expect_equal(unlist(open[1, ]), ls_E_o$o[1], ignore_attr = TRUE)
    expect_equal(unlist(open[2, ]), ls_E_o$inv_o[1], ignore_attr = TRUE)
  }
})


test_that("get_openness returns NA when neither method is usable", {
  lambda <- 50000; chi <- 50000; psi <- 2; alpha <- 10
  ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
  open <- get_openness(lambda, chi, psi, alpha)
  expect_true(is.na(ls_E_o$o[1]))
  expect_true(is.na(ls_E_o$inv_o[1]))
  expect_true(is.na(open[1, ]))
  expect_true(is.na(open[2, ]))
})


test_that("get_openness handles vector input", {
  lambda <- c(-5, 1000, 10000, -8, 50000)
  chi    <- c(10, 2000, 20000, 0.01, 50000)
  psi    <- c(2, 2, 2, 2, 2)
  alpha  <- c(10, 10, 10, 10, 10)

  # Entries where the closed form is usable (small |lambda|).
  ok <- c(TRUE, FALSE, FALSE, TRUE, FALSE)
  ref_o <- mapply(gig_ref, lambda, chi, psi)
  ref_inv_o <- mapply(gig_ref, lambda, chi, psi, MoreArgs = list(inv = TRUE))

  open <- get_openness(lambda, chi, psi, alpha)
  expect_equal(unlist(open[1, ])[ok], ref_o[ok], ignore_attr = TRUE)
  expect_equal(unlist(open[2, ])[ok], ref_inv_o[ok], ignore_attr = TRUE)
})

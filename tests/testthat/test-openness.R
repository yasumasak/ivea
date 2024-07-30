test_that("E_openness", {
  # lambda_enh <- v_alpha_o_e - v_alpha_enh + ls_x$enh_count
  # chi_enh <- 2 * v_alpha_enh * z_enh_actv / ls_x$enh_rlen
  # psi_enh <- 2 * (v_beta_o_e + ls_x$enh_rlen)

  # Normal case
  lambda <- -5
  chi <- 10
  psi <- 2
  alpha <- 10
  # by ghyp
  o_Egig <- ghyp::Egig(lambda, chi, psi, func='x')
  invo_Egig <- ghyp::Egig(lambda, chi, psi, func='1/x')
  # E_openness
  ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
  # get_openness
  expect_open <- get_openness(lambda, chi, psi, alpha)
  # TEST
  expect_equal(ls_E_o$o[1], o_Egig, ignore_attr=T)
  expect_equal(ls_E_o$inv_o[1], invo_Egig, ignore_attr=T)
  expect_equal(unlist(expect_open[1, ]), o_Egig, ignore_attr=T)
  expect_equal(unlist(expect_open[2, ]), invo_Egig, ignore_attr=T)

  # Case that can not be calculated by ghyp package
  lambda <- 1000
  chi <- 2000
  psi <- 2
  alpha <- 10
  # by ghyp
  o_Egig <- ghyp::Egig(lambda, chi, psi, func='x')
  invo_Egig <- ghyp::Egig(lambda, chi, psi, func='1/x')
  # E_openness
  ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
  # get_openness
  expect_open <- get_openness(lambda, chi, psi, alpha)
  # TEST
  expect_true(is.nan(o_Egig))
  expect_true(is.nan(invo_Egig))
  expect_equal(unlist(expect_open[1, ]), ls_E_o$o[1], ignore_attr=T)
  expect_equal(unlist(expect_open[2, ]), ls_E_o$inv_o[1], ignore_attr=T)

  # Case that can not be calculated by ghyp package
  lambda <- 10000
  chi <- 20000
  psi <- 2
  alpha <- 10
  # by ghyp
  o_Egig <- ghyp::Egig(lambda, chi, psi, func='x')
  invo_Egig <- ghyp::Egig(lambda, chi, psi, func='1/x')
  # E_openness
  ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
  # get_openness
  expect_open <- get_openness(lambda, chi, psi, alpha)
  # TEST
  expect_true(is.nan(o_Egig))
  expect_true(is.nan(invo_Egig))
  expect_equal(unlist(expect_open[1, ]), ls_E_o$o[1], ignore_attr=T)
  expect_equal(unlist(expect_open[2, ]), ls_E_o$inv_o[1], ignore_attr=T)

  # Case that can not be calculated by E_openness
  lambda <- -8
  chi <- 0.01
  psi <- 2
  alpha <- 10
  # by ghyp
  o_Egig <- ghyp::Egig(lambda, chi, psi, func='x')
  invo_Egig <- ghyp::Egig(lambda, chi, psi, func='1/x')
  # E_openness
  ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
  # get_openness
  expect_open <- get_openness(lambda, chi, psi, alpha)
  # TEST
  expect_true(is.na(ls_E_o$o[1]))
  expect_true(is.na(ls_E_o$inv_o[1]))
  expect_equal(unlist(expect_open[1, ]), o_Egig, ignore_attr=T)
  expect_equal(unlist(expect_open[2, ]), invo_Egig, ignore_attr=T)

  # Case that can not be calculated
  lambda <- 50000
  chi <- 50000
  psi <- 2
  alpha <- 10
  # by ghyp
  o_Egig <- ghyp::Egig(lambda, chi, psi, func='x')
  invo_Egig <- ghyp::Egig(lambda, chi, psi, func='1/x')
  # E_openness
  ls_E_o <- IVEA::E_openness(lambda, chi, psi, alpha)
  # get_openness
  expect_open <- get_openness(lambda, chi, psi, alpha)
  # TEST
  expect_true(is.nan(o_Egig))
  expect_true(is.nan(invo_Egig))
  expect_true(is.na(ls_E_o$o[1]))
  expect_true(is.na(ls_E_o$inv_o[1]))
  expect_true(is.na(expect_open[1, ]))
  expect_true(is.na(expect_open[2, ]))

})



test_that("E_openness vector input", {
  # lambda_enh <- v_alpha_o_e - v_alpha_enh + ls_x$enh_count
  # chi_enh <- 2 * v_alpha_enh * z_enh_actv / ls_x$enh_rlen
  # psi_enh <- 2 * (v_beta_o_e + ls_x$enh_rlen)

  # Parameters in vectors
  lambda <- c(-5, 1000, 10000, -8, 50000)
  chi <- c(10, 2000, 20000, 0.01, 50000)
  psi <- c(2, 2, 2, 2, 2)
  alpha <- c(10, 10, 10, 10, 10)

  # Calculated by ghyp
  v_ok_ghyp <- c(T, F, F, T, F)

  # by ghyp
  o_Egig <- ghyp::Egig(lambda, chi, psi, func='x')
  invo_Egig <- ghyp::Egig(lambda, chi, psi, func='1/x')
  # get_openness
  expect_open <- get_openness(lambda, chi, psi, alpha)

  # TEST
  expect_equal(o_Egig[v_ok_ghyp], unlist(expect_open[1, ])[v_ok_ghyp], ignore_attr=T)
  expect_equal(invo_Egig[v_ok_ghyp], unlist(expect_open[2, ])[v_ok_ghyp], ignore_attr=T)

})

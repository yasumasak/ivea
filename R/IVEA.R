
#' Run variational Bayesian inference
#'
#' Perform variational inference of regulatory elements (promoter and enhancers)
#' opennesses and activities.
#'
#' @details
#' The parameters and expectations of the latent variables are iteratively updated
#' using the equations described in (reference) until convergence of the enhancer
#' activities judged by `stop_criterion` or `max_iteration`.
#'
#' @param param a parameter list object containing the following components:
#'   \describe{
#'     \item{alpha_pro}{Hyperparameter alpha for the gamma prior of promoter activity.}
#'     \item{alpha_enh}{Hyperparameter alpha for the gamma prior of enhancer activity.}
#'     \item{alpha_o}{Hyperparameter alpha for the gamma prior of openness.}
#'     \item{beta_o}{Hyperparameter beta for the gamma prior of openness.}
#'     \item{alpha_k}{Hyperparameter alpha for the gamma prior of skaler k.}
#'     \item{beta_k}{Hyperparameter beta for the gamma prior of skaler k.}
#'     \item{no_le_enh}{use IVEA_nolE model or not.}
#'     \item{max_iteration}{Maximum iterations in variational inference.}
#'     \item{stop_criterion}{Stop iteration when all relative change of
#'                           enhancer activities are less than this value.}
#'   }
#' @param ls_x an input-data list object containing the following components:
#'   \describe{
#'     \item{chr}{chromosome.}
#'     \item{region}{region on chromosome.}
#'     \item{gene_name}{character vector of gene names.}
#'     \item{gene_tss}{numeric vector of gene TSSs.}
#'     \item{rna_count}{numeric vector of sequence read counts on genes.}
#'     \item{rna_rlen}{numeric vector of relative effective lengths of genes.}
#'     \item{pro_name}{character vector of gene promoter names.}
#'     \item{pro_count}{numeric vecor of sequence read counts gene promoters.}
#'     \item{pro_rlen}{numeric vecor of relative lengths of gene promoters.}
#'     \item{enh_name}{character vector of enhancer names.}
#'     \item{enh_count}{numeric vecor of sequence read counts enhancers.}
#'     \item{enh_rlen}{numeric vecor of relative lengths of enhancers.}
#'     \item{n_gene}{total number of genes.}
#'     \item{n_enh}{total number of enhancers.}
#'     \item{b_size}{numeric vector of genomic sequence-based burst size.}
#'     \item{contact}{matrix (genes x enhancers) of contact frequencies.}
#'     \item{sp_contact}{sparse matrix (genes x enhancers) of contact frequencies.}
#'     \item{sp_distance}{sparse matrix (genes x enhancers) of genomic distances.}
#'   }
#'
#' @return A resultant list object containing the following components:
#'   \describe{
#'     \item{k}{vector of estimated scaling factor values.}
#'     \item{pro_openness}{matrix of estimated opennesses of gene promoters.}
#'     \item{enh_openness}{matrix of estimated opennesses of enhancers.}
#'     \item{pro_inv_openness}{matrix of estimated 1/opennesses of gene promoters.}
#'     \item{enh_inv_openness}{matrix of estimated 1/opennesses of enhancers.}
#'     \item{pro_activity}{matrix of estimated activities of gene promoters.}
#'     \item{enh_activity}{matrix of estimated activities of enhancers.}
#'     \item{pro_log_activity}{matrix of estimated log activities of gene promoters.}
#'     \item{enh_log_activity}{matrix of estimated log activities of enhancers.}
#'   }
#'
#' @export
run_variational_bayes <- function(param, ls_x){

  #------------------------------
  # Matrices to store values from the iterations
  n_iter <- param$max_iteration + 1

  #------------------------------
  # Hyperprameters
  v_alpha_pro <- rep(param$alpha_pro, ls_x$n_gene)
  v_alpha_enh <- rep(param$alpha_enh, ls_x$n_enh)
  v_alpha_o_p <- rep(param$alpha_o, ls_x$n_gene)
  v_alpha_o_e <- rep(param$alpha_o, ls_x$n_enh)
  v_beta_o_p <- rep(param$beta_o, ls_x$n_gene)
  v_beta_o_e <- rep(param$beta_o, ls_x$n_enh)
  alpha_k <- param$alpha_k
  beta_k <- param$beta_k

  #------------------------------
  # Make a resultant object with initial values
  ls_z <- make_initial_z(ls_x, n_iter, param$no_le_enh)

  #------------------------------
  # Valuables to be updated
  z_k <- ls_z$k[1]
  z_pro_open <- ls_z$pro_openness[1, ]
  z_enh_open <- ls_z$enh_openness[1, ]
  z_pro_actv <- ls_z$pro_activity[1, ]
  z_enh_actv <- ls_z$enh_activity[1, ]
  z_pro_inv_open <- ls_z$pro_inv_openness[1, ]
  z_enh_inv_open <- ls_z$enh_inv_openness[1, ]
  z_pro_log_actv <- ls_z$pro_log_activity[1, ]
  z_enh_log_actv <- ls_z$enh_log_activity[1, ]

  # For checking convergence
  is_converged <- F
  v_max_rl_change <- vector(length=n_iter)

  cat("  Start iteration... ")
  #------------------------------
  # Iteration of variational Bayes
  z_k <- ls_z$k[1]
  for (iter in 2:n_iter) {

    # Update promoter openness ------------------------------
    # Generalized inverse Gaussian distribution parameters
    lambda_pro <- v_alpha_o_p - v_alpha_pro + ls_x$pro_count
    chi_pro <- 2 * v_alpha_pro * z_pro_actv / ls_x$b_size
    psi_pro <- 2 * (v_beta_o_p + ls_x$pro_rlen)
    # Expected values
    expect_open_pro <- get_openness(lambda=lambda_pro, chi=chi_pro, psi=psi_pro, alpha=v_alpha_pro)
    o_pro <- unlist(expect_open_pro[1, ])
    inv_o_pro <- unlist(expect_open_pro[2, ])
    z_pro_open[!is.na(o_pro)] <- o_pro[!is.na(o_pro)]
    z_pro_inv_open[!is.na(inv_o_pro)] <- inv_o_pro[!is.na(inv_o_pro)]

    # Update enhancer openness ------------------------------
    # Generalized inverse Gaussian distribution parameters
    lambda_enh <- v_alpha_o_e - v_alpha_enh + ls_x$enh_count
    if(param$no_le_enh){
      # IVEA_nolE
      chi_enh <- 2 * v_alpha_enh * z_enh_actv / 1
    }else{
      # IVEA
      chi_enh <- 2 * v_alpha_enh * z_enh_actv / ls_x$enh_rlen
    }
    psi_enh <- 2 * (v_beta_o_e + ls_x$enh_rlen)
    # Expected values
    expect_open_enh <- get_openness(lambda=lambda_enh, chi=chi_enh, psi=psi_enh, alpha=v_alpha_enh)
    o_enh <- unlist(expect_open_enh[1, ])
    inv_o_enh <- unlist(expect_open_enh[2, ])
    z_enh_open[!is.na(o_enh)] <- o_enh[!is.na(o_enh)]
    z_enh_inv_open[!is.na(inv_o_enh)] <- inv_o_enh[!is.na(inv_o_enh)]

    # Update promoter activity ------------------------------
    # Gamma distribution parameters
    a_pro_actv <- v_alpha_pro + ls_x$rna_count
    b_pro_actv <- (v_alpha_pro * z_pro_inv_open / ls_x$b_size
                   + drop(ls_x$rna_rlen * z_k * (ls_x$contact %*% z_enh_actv)))
    # Expected values
    z_pro_actv <- a_pro_actv / b_pro_actv
    z_pro_log_actv <- digamma(a_pro_actv) - log(b_pro_actv)

    # Update enhancer activity ------------------------------
    # Auxiliary probability
    v_sum_aux <- ls_x$contact %*% exp(z_enh_log_actv)
    mx_aux <- sweep(sweep(ls_x$contact, 2, exp(z_enh_log_actv), '*'), 1, v_sum_aux, '/')
    # Gamma distribution parameters
    a_enh_actv <- v_alpha_enh + drop(ls_x$rna_count %*% mx_aux)
    if(param$no_le_enh){
      # IVEA_nolE
      b_enh_actv <- (v_alpha_enh * z_enh_inv_open / 1
                     + colSums( sweep(ls_x$contact, 1, (ls_x$rna_rlen * z_k * z_pro_actv), '*') ))
    }else{
      # IVEA
      b_enh_actv <- (v_alpha_enh * z_enh_inv_open / ls_x$enh_rlen
                     + colSums( sweep(ls_x$contact, 1, (ls_x$rna_rlen * z_k * z_pro_actv), '*') ))
    }
    # Expected values
    z_enh_actv <- a_enh_actv / b_enh_actv
    z_enh_log_actv <- digamma(a_enh_actv) - log(b_enh_actv)

    # Update Scaling factor k ------------------------------
    # Gamma distribution parameters
    a_k <- alpha_k + sum(ls_x$rna_count)
    b_k <- beta_k + sum(ls_x$rna_rlen * z_pro_actv * (ls_x$contact %*% z_enh_actv))
    # Expected value
    z_k <- a_k / b_k

    # Store updated values ------------------------------
    ls_z$k[iter] <- z_k
    ls_z$pro_openness[iter, ] <- z_pro_open
    ls_z$enh_openness[iter, ] <- z_enh_open
    ls_z$pro_activity[iter, ] <- z_pro_actv
    ls_z$pro_activity_shape[iter, ] <- a_pro_actv
    ls_z$pro_activity_rate[iter, ] <- b_pro_actv
    ls_z$enh_activity[iter, ] <- z_enh_actv
    ls_z$enh_activity_shape[iter, ] <- a_enh_actv
    ls_z$enh_activity_rate[iter, ] <- b_enh_actv
    ls_z$pro_inv_openness[iter, ] <- z_pro_inv_open
    ls_z$enh_inv_openness[iter, ] <- z_enh_inv_open
    ls_z$pro_log_activity[iter, ] <- z_pro_log_actv
    ls_z$enh_log_activity[iter, ] <- z_enh_log_actv

    # Check convergence ------------------------------
    v_rl_change <- (abs(ls_z$enh_activity[iter, ] - ls_z$enh_activity[iter - 1, ])
                    / ls_z$enh_activity[iter, ])
    max_rl_change <- max(v_rl_change)
    v_max_rl_change[iter] <- max_rl_change
    if(!is.na(max_rl_change) & max_rl_change < param$stop_criterion){
      is_converged <- T
    }

    # When converged ------------------------------
    if(is_converged | iter==n_iter){
      if(is_converged){
        cat("Converged at ",iter - 1,"\n")
      }else{
        cat("Finished at max iteration ",iter - 1,"\n")
      }
      break
    }

    if(iter %% 10 == 0){ cat(iter," ") } # Keep track of progress
  }

  ls_z$k <- ls_z$k[1:iter]
  ls_z$pro_openness <- ls_z$pro_openness[1:iter, ]
  ls_z$enh_openness <- ls_z$enh_openness[1:iter, ]
  ls_z$pro_activity <- ls_z$pro_activity[1:iter, ]
  ls_z$pro_activity_shape <- ls_z$pro_activity_shape[1:iter, ]
  ls_z$pro_activity_rate <- ls_z$pro_activity_rate[1:iter, ]
  ls_z$enh_activity <- ls_z$enh_activity[1:iter, ]
  ls_z$enh_activity_shape <- ls_z$enh_activity_shape[1:iter, ]
  ls_z$enh_activity_rate <- ls_z$enh_activity_rate[1:iter, ]
  ls_z$pro_inv_openness <- ls_z$pro_inv_openness[1:iter, ]
  ls_z$enh_inv_openness <- ls_z$enh_inv_openness[1:iter, ]
  ls_z$pro_log_activity <- ls_z$pro_log_activity[1:iter, ]
  ls_z$enh_log_activity <- ls_z$enh_log_activity[1:iter, ]
  return(ls_z)

}


#' Make an initial resultant object
#'
#' Make a list object for latent variables to be estimated and set the initial values
#'
#' @param ls_x an input-data list object containing the following components:
#'   \describe{
#'     \item{chr}{chromosome.}
#'     \item{region}{region on chromosome.}
#'     \item{gene_name}{character vector of gene names.}
#'     \item{gene_tss}{numeric vector of gene TSSs.}
#'     \item{rna_count}{numeric vector of sequence read counts on genes.}
#'     \item{rna_rlen}{numeric vector of relative effective lengths of genes.}
#'     \item{pro_name}{character vector of gene promoter names.}
#'     \item{pro_count}{numeric vecor of sequence read counts gene promoters.}
#'     \item{pro_rlen}{numeric vecor of relative lengths of gene promoters.}
#'     \item{enh_name}{character vector of enhancer names.}
#'     \item{enh_count}{numeric vecor of sequence read counts enhancers.}
#'     \item{enh_rlen}{numeric vecor of relative lengths of enhancers.}
#'     \item{n_gene}{total number of genes.}
#'     \item{n_enh}{total number of enhancers.}
#'     \item{b_size}{numeric vector of genomic sequence-based burst size.}
#'     \item{contact}{matrix (genes x enhancers) of contact frequencies.}
#'     \item{sp_contact}{sparse matrix (genes x enhancers) of contact frequencies.}
#'     \item{sp_distance}{sparse matrix (genes x enhancers) of genomic distances.}
#'   }
#' @param n_iter Maximum number of iterations
#' @param no_le_enh Use IVEA_nolE model or not
#'
#' @return An initial resultant list object containing the following components:
#'   \describe{
#'     \item{k}{vector of estimated scaling factor values.}
#'     \item{pro_openness}{matrix of estimated opennesses of gene promoters.}
#'     \item{enh_openness}{matrix of estimated opennesses of enhancers.}
#'     \item{pro_inv_openness}{matrix of estimated 1/opennesses of gene promoters.}
#'     \item{enh_inv_openness}{matrix of estimated 1/opennesses of enhancers.}
#'     \item{pro_activity}{matrix of estimated activities of gene promoters.}
#'     \item{enh_activity}{matrix of estimated activities of enhancers.}
#'     \item{pro_log_activity}{matrix of estimated log activities of gene promoters.}
#'     \item{enh_log_activity}{matrix of estimated log activities of enhancers.}
#'   }
make_initial_z <- function(ls_x, n_iter, no_le_enh){

  # Vector and matrix to store estimated values
  v_k <- vector(length=n_iter)
  pro_openness <- matrix(nrow=n_iter, ncol=ls_x$n_gene)
  enh_openness <- matrix(nrow=n_iter, ncol=ls_x$n_enh)
  pro_inv_openness <- matrix(nrow=n_iter, ncol=ls_x$n_gene)
  enh_inv_openness <- matrix(nrow=n_iter, ncol=ls_x$n_enh)
  pro_activity <- matrix(nrow=n_iter, ncol=ls_x$n_gene)
  pro_activity_shape <- matrix(nrow=n_iter, ncol=ls_x$n_gene)
  pro_activity_rate <- matrix(nrow=n_iter, ncol=ls_x$n_gene)
  enh_activity <- matrix(nrow=n_iter, ncol=ls_x$n_enh)
  enh_activity_shape <- matrix(nrow=n_iter, ncol=ls_x$n_enh)
  enh_activity_rate <- matrix(nrow=n_iter, ncol=ls_x$n_enh)
  pro_log_activity <- matrix(nrow=n_iter, ncol=ls_x$n_gene)
  enh_log_activity <- matrix(nrow=n_iter, ncol=ls_x$n_enh)

  # Column names of the matrices
  colnames(pro_openness) <- ls_x$gene_name
  colnames(enh_openness) <- ls_x$enh_name
  colnames(pro_inv_openness) <- ls_x$gene_name
  colnames(enh_inv_openness) <- ls_x$enh_name
  colnames(pro_activity) <- ls_x$gene_name
  colnames(pro_activity_shape) <- ls_x$gene_name
  colnames(pro_activity_rate) <- ls_x$gene_name
  colnames(enh_activity) <- ls_x$enh_name
  colnames(enh_activity_shape) <- ls_x$enh_name
  colnames(enh_activity_rate) <- ls_x$enh_name
  colnames(pro_log_activity) <- ls_x$gene_name
  colnames(enh_log_activity) <- ls_x$enh_name

  # List of the resultant vector and matrices
  ls_z <- list(k=v_k,
               pro_openness=pro_openness,
               enh_openness=enh_openness,
               pro_inv_openness=pro_inv_openness,
               enh_inv_openness=enh_inv_openness,
               pro_activity=pro_activity,
               pro_activity_shape=pro_activity_shape,
               pro_activity_rate=pro_activity_rate,
               enh_activity=enh_activity,
               enh_activity_shape=enh_activity_shape,
               enh_activity_rate=enh_activity_rate,
               pro_log_activity=pro_log_activity,
               enh_log_activity=enh_log_activity)

  # Set initial values
  ls_z$k[1] <- ( (2/ls_x$n_gene)
                 * sum( (ls_x$rna_count / ls_x$rna_rlen) /
                          ((ls_x$b_size * ls_x$pro_count / ls_x$pro_rlen)
                            * (ls_x$contact %*% ls_x$enh_count)) ) )
  ls_z$pro_openness[1,] <- ls_x$pro_count / ls_x$pro_rlen + 0.01
  ls_z$enh_openness[1,] <- ls_x$enh_count / ls_x$enh_rlen + 0.01
  ls_z$pro_inv_openness[1,] <- ls_x$pro_rlen / (ls_x$pro_count + 0.01)
  ls_z$enh_inv_openness[1,] <- ls_x$enh_rlen / (ls_x$enh_count + 0.01)
  ls_z$pro_activity[1, ] <- ls_x$pro_count * ls_x$b_size / ls_x$pro_rlen + 0.01
  if(no_le_enh){
    # IVEA_nolE
    ls_z$enh_activity[1, ] <- ls_x$enh_count / ls_x$enh_rlen + 0.01
  }else{
    # IVEA
    ls_z$enh_activity[1, ] <- ls_x$enh_count + 0.01
  }
  ls_z$pro_log_activity[1,] <- log(ls_x$pro_count * ls_x$b_size / ls_x$pro_rlen + 0.01)
  ls_z$enh_log_activity[1,] <- log(ls_x$enh_count + 0.01)

  return(ls_z)
}


#' Make summary data frame
#'
#' Calculate enhancer-gene pairs' strength, contribution and score,
#' and make summary data frames of enhancer-gene pairs
#'
#' @param ls_x an input-data list object containing the following components:
#'   \describe{
#'     \item{chr}{chromosome.}
#'     \item{region}{region on chromosome.}
#'     \item{gene_name}{character vector of gene names.}
#'     \item{gene_tss}{numeric vector of gene TSSs.}
#'     \item{rna_count}{numeric vector of sequence read counts on genes.}
#'     \item{rna_rlen}{numeric vector of relative effective lengths of genes.}
#'     \item{pro_name}{character vector of gene promoter names.}
#'     \item{pro_count}{numeric vecor of sequence read counts gene promoters.}
#'     \item{pro_rlen}{numeric vecor of relative lengths of gene promoters.}
#'     \item{enh_name}{character vector of enhancer names.}
#'     \item{enh_count}{numeric vecor of sequence read counts enhancers.}
#'     \item{enh_rlen}{numeric vecor of relative lengths of enhancers.}
#'     \item{n_gene}{total number of genes.}
#'     \item{n_enh}{total number of enhancers.}
#'     \item{b_size}{numeric vector of genomic sequence-based burst size.}
#'     \item{contact}{matrix (genes x enhancers) of contact frequencies.}
#'     \item{sp_contact}{sparse matrix (genes x enhancers) of contact frequencies.}
#'     \item{sp_distance}{sparse matrix (genes x enhancers) of genomic distances.}
#'   }
#' @param ls_z a resultant list object containing the following components:
#'   \describe{
#'     \item{k}{vector of estimated scaling factor values.}
#'     \item{pro_openness}{matrix of estimated opennesses of gene promoters.}
#'     \item{enh_openness}{matrix of estimated opennesses of enhancers.}
#'     \item{pro_inv_openness}{matrix of estimated 1/opennesses of gene promoters.}
#'     \item{enh_inv_openness}{matrix of estimated 1/opennesses of enhancers.}
#'     \item{pro_activity}{matrix of estimated activities of gene promoters.}
#'     \item{enh_activity}{matrix of estimated activities of enhancers.}
#'     \item{pro_log_activity}{matrix of estimated log activities of gene promoters.}
#'     \item{enh_log_activity}{matrix of estimated log activities of enhancers.}
#'   }
#'
#' @return A list object containing the following summary data frames:
#'   \describe{
#'     \item{info_pair}{data frame of information of enhancer-gene pairs.}
#'     \item{bedpe_pair}{data frame of prediction scores in BEDPE format.}
#'   }
#'
#' @export
make_summary <- function(ls_x, ls_z){

  # Estimated element activities
  z_pro_actv <- as.vector(utils::tail(ls_z$pro_activity, n=1))
  z_pro_actv_shape <- as.vector(utils::tail(ls_z$pro_activity_shape, n=1))
  z_pro_actv_rate <- as.vector(utils::tail(ls_z$pro_activity_rate, n=1))
  z_enh_actv <- as.vector(utils::tail(ls_z$enh_activity, n=1))
  z_enh_actv_shape <- as.vector(utils::tail(ls_z$enh_activity_shape, n=1))
  z_enh_actv_rate <- as.vector(utils::tail(ls_z$enh_activity_rate, n=1))

  # Strength / Contribution / Score of enhancer-gene pairs
  mx_contr <- ls_x$contact * (z_pro_actv %*% t(z_enh_actv))
  mx_stren <- ls_x$contact * (rep(1, ls_x$n_gene) %*% t(z_enh_actv))
  mx_score <- mx_stren / (rowSums(mx_stren) %*% t(rep(1, ls_x$n_enh)))

  # Summary of sparse matrix
  ssp_stren <- Matrix::summary(methods::as(mx_stren, 'sparseMatrix'))
  ssp_contr <- Matrix::summary(methods::as(mx_contr, 'sparseMatrix'))
  ssp_score <- Matrix::summary(methods::as(mx_score, 'sparseMatrix'))
  ssp_dist <- Matrix::summary(ls_x$sp_distance)
  ssp_con <- Matrix::summary(ls_x$sp_contact)

  if(sum(ssp_con[, c(1,2)] != ssp_score[, c(1,2)]) > 0){stop("indeces are different")}

  # Pairs' information
  n_pairs <- nrow(ssp_score)
  df_info_pair <- data.frame(gene=ls_x$gene_name[ssp_score$i],
                             gene_chr=rep(ls_x$chr, n_pairs),
                             gene_tss=ls_x$gene_tss[ssp_score$i],
                             promoter=ls_x$pro_name[ssp_score$i],
                             promoter_count=ls_x$pro_count[ssp_score$i],
                             promoter_length=ls_x$pro_rlen[ssp_score$i],
                             promoter_activity=z_pro_actv[ssp_score$i],
                             enhancer=ls_x$enh_name[ssp_score$j],
                             enhancer_count=ls_x$enh_count[ssp_score$j],
                             enhancer_length=ls_x$enh_rlen[ssp_score$j],
                             enhancer_activity=z_enh_actv[ssp_score$j],
                             distance=ssp_dist$x,
                             contact=ssp_con$x,
                             strength=ssp_stren$x,
                             contribution=ssp_contr$x,
                             score=ssp_score$x)

  # Pairs' score in BEDPE format
  v_pair_name <- paste(ls_x$gene_name[ssp_score$i], ls_x$enh_name[ssp_score$j], sep="-")
  mx_e <- matrix(unlist(strsplit(ls_x$enh_name[ssp_score$j], split="[:-]")), ncol=3, byrow=T)
  df_bedpe_pair <- data.frame(e_chr=mx_e[, 1],
                              e_start=mx_e[, 2],
                              e_end=mx_e[, 3],
                              g_chr=rep(ls_x$chr, n_pairs),
                              g_start=ls_x$gene_tss[ssp_score$i]-1,
                              g_end=ls_x$gene_tss[ssp_score$i],
                              name=v_pair_name, score=ssp_score$x,
                              e_strand=rep(".", n_pairs),
                              g_strand=rep(".", n_pairs))


  # Promoter activity in in BED format
  z_pro_actv_025 <- qgamma(0.025, shape = z_pro_actv_shape, rate = z_pro_actv_rate)
  z_pro_actv_050 <- qgamma(0.5, shape = z_pro_actv_shape, rate = z_pro_actv_rate)
  z_pro_actv_975 <- qgamma(0.975, shape = z_pro_actv_shape, rate = z_pro_actv_rate)
  v_sd_pro_actv <- sqrt(z_pro_actv_shape/(z_pro_actv_rate^2))
  df_bed_pro <- data.frame(p_chr=ls_x$chr,
                           p_start=ls_x$gene_tss-1,
                           p_end=ls_x$gene_tss,
                           name=ls_x$gene_name,
                           activity=z_pro_actv,
                           activity_025=z_pro_actv_025,
                           activity_050=z_pro_actv_050,
                           activity_975=z_pro_actv_975,
                           sd_activity=v_sd_pro_actv,
                           shape=z_pro_actv_shape,
                           rate=z_pro_actv_rate
                           )

  # Enhancer activity in BED format
  mx_enh <- matrix(unlist(strsplit(ls_x$enh_name, split="[:-]")), ncol=3, byrow=T)
  z_enh_actv_025 <- qgamma(0.025, shape = z_enh_actv_shape, rate = z_enh_actv_rate)
  z_enh_actv_050 <- qgamma(0.5, shape = z_enh_actv_shape, rate = z_enh_actv_rate)
  z_enh_actv_975 <- qgamma(0.975, shape = z_enh_actv_shape, rate = z_enh_actv_rate)
  v_sd_enh_actv <- sqrt(z_enh_actv_shape/(z_enh_actv_rate^2))
  df_bed_enh <- data.frame(e_chr=mx_enh[, 1],
                           e_start=mx_enh[, 2],
                           e_end=mx_enh[, 3],
                           name=ls_x$enh_name,
                           activity=z_enh_actv,
                           activity_025=z_enh_actv_025,
                           activity_050=z_enh_actv_050,
                           activity_975=z_enh_actv_975,
                           sd_activity=v_sd_enh_actv,
                           shape=z_enh_actv_shape,
                           rate=z_enh_actv_rate
                           )

  # Order by score
  order_score <- order(df_bedpe_pair$score, decreasing=T)
  df_info_pair <- df_info_pair[order_score, ]
  df_bedpe_pair <- df_bedpe_pair[order_score, ]

  return(list(info_pair=df_info_pair, bedpe_pair=df_bedpe_pair,
               bed_pro=df_bed_pro, bed_enh=df_bed_enh))

}


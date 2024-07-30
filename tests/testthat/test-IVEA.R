get_opts <- function(){
  opts <- list(
    chr_regions= "chr22",
    alpha_pro = 80,
    alpha_enh = 10,
    alpha_o = 0.01,
    beta_o = 0.0001,
    alpha_k = 0.001,
    beta_k = 0.001,
    # Iteration
    max_iteration = 500,
    stop_criterion = 0.005,
    # Chromatin contact
    contacts = "../../example/output/enhancer-gene_contacts.chr22.txt",
    hic_contact_colname = NULL,
    window = 3000000,
    # Gene expression
    rnas = "../../example/output/gene_expressions.txt",
    rna_id_colname = "gene_name",
    rna_length_colname = "effective_length",
    rna_count_colname = "expected_count",
    rna_cutoff_colname = "TPM",
    rna_cutoff_value = 8,
    rna_type = "gene",
    # Regulatory element
    enhancers = "../../example/output/enhancer_elements.txt",
    promoters = "../../example/output/gene_promoter_elements.txt",
    burst_sizes = "../../reference/gene_burst_sizes.hg19.txt"
  )
  return(opts)
}

test_that("IVEA example chr22", {
  # Options
  opts <- get_opts()
  # Original IVEA
  opts$no_le_enh <- FALSE

  chr <- opts$chr_regions
  region <- NULL
  label <- opts$chr_regions
  rlen_unit=1000

  # load_data
  ls_x <- IVEA::load_data(opts, chr, region)

  # Tested function: run_variational_bayes
  ls_z <- IVEA::run_variational_bayes(opts, ls_x)
  n_iter <- length(ls_z$k)
  # TEST
  expect_equal(dim(ls_z$pro_openness), c(n_iter, ls_x$n_gene), ignore_attr=T)
  expect_equal(dim(ls_z$pro_inv_openness), c(n_iter, ls_x$n_gene), ignore_attr=T)
  expect_equal(dim(ls_z$pro_activity), c(n_iter, ls_x$n_gene), ignore_attr=T)
  expect_equal(dim(ls_z$pro_log_activity), c(n_iter, ls_x$n_gene), ignore_attr=T)
  expect_equal(dim(ls_z$enh_openness), c(n_iter, ls_x$n_enh), ignore_attr=T)
  expect_equal(dim(ls_z$enh_inv_openness), c(n_iter, ls_x$n_enh), ignore_attr=T)
  expect_equal(dim(ls_z$enh_activity), c(n_iter, ls_x$n_enh), ignore_attr=T)
  expect_equal(dim(ls_z$enh_log_activity), c(n_iter, ls_x$n_enh), ignore_attr=T)
  expect_equal(sum(is.na(ls_z$k)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$pro_openness)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$pro_inv_openness)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$pro_activity)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$pro_log_activity)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$enh_openness)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$enh_inv_openness)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$enh_activity)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$enh_log_activity)), 0, ignore_attr=T)

  # Tested function: make_summary
  ls_out <- IVEA::make_summary(ls_x, ls_z)

  # Check Pair information with valid results
  df_bedpe_pair <- read.table("../../example/output/predictions_score.chr22.bedpe", sep="\t", header=F)
  df_info_pair <- read.table("../../example/output/predictions_info.chr22.txt", sep="\t", header=T)

  # Pair name
  expect_equal(ls_out$bedpe_pair[, 7], df_bedpe_pair[, 7], ignore_attr=T)
  # Score
  expect_equal(signif(ls_out$bedpe_pair[, 8], 6), signif(df_bedpe_pair[, 8], 6), ignore_attr=T)
  # Promoter/Enhancer activity
  expect_equal(signif(ls_out$info_pair[, "promoter_activity"], 6), signif(df_info_pair[, "promoter_activity"], 6), ignore_attr=T)
  expect_equal(signif(ls_out$info_pair[, "enhancer_activity"], 6), signif(df_info_pair[, "enhancer_activity"], 6), ignore_attr=T)
  # Score
  expect_equal(signif(ls_out$info_pair[, "score"], 6), signif(df_info_pair[, "score"], 6), ignore_attr=T)
  expect_equal(signif(ls_out$bedpe_pair[, 8], 6), signif(ls_out$info_pair[, "score"], 6), ignore_attr=T)

  # Check estimated Promoter/Enhancer activity with valid results
  df_pro_actv <- read.table("../../example/output/estimates.promoter_activity.chr22.bed", sep="\t", header=F)
  df_enh_actv <- read.table("../../example/output/estimates.enhancer_activity.chr22.bed", sep="\t", header=F)

  # TEST
  expect_equal(nrow(ls_out$bed_pro), ls_x$n_gene, ignore_attr=T)
  expect_equal(nrow(ls_out$bed_enh), ls_x$n_enh, ignore_attr=T)
  expect_equal(ls_out$bed_pro[, 4], df_pro_actv[, 4], ignore_attr=T)
  expect_equal(ls_out$bed_enh[, 4], df_enh_actv[, 4], ignore_attr=T)
  expect_equal(signif(ls_out$bed_pro[, 5], 6), signif(df_pro_actv[, 5], 6), ignore_attr=T)
  expect_equal(signif(ls_out$bed_enh[, 5], 6), signif(df_enh_actv[, 5], 6), ignore_attr=T)

  # Confirm output promoter activities are the same
  df_pro_bed <- unique(ls_out$info_pair[, c("gene", "promoter_activity")])
  rownames(df_pro_bed) <- df_pro_bed$gene
  expect_equal(df_pro_bed[ls_x$gene_name, ], ls_out$bed_pro[, c("name", "activity")], ignore_attr=T)
  # Confirm output enhancer activities are the same
  df_enh_bed <- unique(ls_out$info_pair[, c("enhancer", "enhancer_activity")])
  rownames(df_enh_bed) <- df_enh_bed$enhancer
  expect_equal(df_enh_bed[ls_x$enh_name, ], ls_out$bed_enh[, c("name", "activity")], ignore_attr=T)
})



test_that("IVEA_nolE example chr22", {
  # Options
  opts <- get_opts()
  # IVEA_nolE
  # omitting the relative enhancer lengths in the prior information of enhancer activities
  opts$no_le_enh <- TRUE

  chr <- opts$chr_regions
  region <- NULL
  label <- opts$chr_regions
  rlen_unit=1000

  # load_data
  ls_x <- IVEA::load_data(opts, chr, region)

  # Tested function: run_variational_bayes
  ls_z <- IVEA::run_variational_bayes(opts, ls_x)
  ls_out <- IVEA::make_summary(ls_x, ls_z)
  n_iter <- length(ls_z$k)

  # TEST
  expect_equal(dim(ls_z$pro_activity), c(n_iter, ls_x$n_gene), ignore_attr=T)
  expect_equal(dim(ls_z$enh_activity), c(n_iter, ls_x$n_enh), ignore_attr=T)
  expect_equal(sum(is.na(ls_z$pro_activity)), 0, ignore_attr=T)
  expect_equal(sum(is.na(ls_z$enh_activity)), 0, ignore_attr=T)

  # Check Pair information with valid results
  df_bedpe_pair <- read.table("../../example/output/predictions_score.chr22.bedpe", sep="\t", header=F)
  df_info_pair <- read.table("../../example/output/predictions_info.chr22.txt", sep="\t", header=T)

  # Score (should be different from the original IVEA result)
  expect_false(any(signif(ls_out$bedpe_pair[, 8], 6) == signif(df_bedpe_pair[, 8], 6)))
  # Promoter/Enhancer activity (should be different from the original IVEA result)
  expect_false(any(signif(ls_out$info_pair[, "promoter_activity"], 6) == signif(df_info_pair[, "promoter_activity"], 6)))
  expect_false(any(signif(ls_out$info_pair[, "enhancer_activity"], 6) == signif(df_info_pair[, "enhancer_activity"], 6)))

})


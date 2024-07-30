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

test_that("data example chr22", {
  # Options
  opts <- get_opts()
  chr <- opts$chr_regions
  region <- NULL
  label <- opts$chr_regions
  rlen_unit=1000

  # Tested function: load_data
  ls_x <- IVEA::load_data(opts, chr, region)

  # TEST
  expect_equal(ls_x$chr, "chr22", ignore_attr=T)
  expect_null(ls_x$region)
  expect_equal(ls_x$n_gene, 193, ignore_attr=T)
  expect_equal(ls_x$n_enh, 3174, ignore_attr=T)
  expect_equal(names(ls_x$b_size), ls_x$gene_name, ignore_attr=T)
  expect_equal(length(ls_x$gene_tss), ls_x$n_gene, ignore_attr=T)
  expect_equal(length(ls_x$rna_count), ls_x$n_gene, ignore_attr=T)
  expect_equal(length(ls_x$rna_rlen), ls_x$n_gene, ignore_attr=T)
  expect_equal(length(ls_x$pro_name), ls_x$n_gene, ignore_attr=T)
  expect_equal(length(ls_x$pro_count), ls_x$n_gene, ignore_attr=T)
  expect_equal(length(ls_x$pro_rlen), ls_x$n_gene, ignore_attr=T)
  expect_equal(length(ls_x$enh_name), ls_x$n_enh, ignore_attr=T)
  expect_equal(length(ls_x$enh_count), ls_x$n_enh, ignore_attr=T)
  expect_equal(length(ls_x$enh_rlen), ls_x$n_enh, ignore_attr=T)
  expect_equal(dim(ls_x$contact), c(ls_x$n_gene, ls_x$n_enh), ignore_attr=T)
  expect_gte(min(ls_x$contact), 0)
  expect_lt(max(ls_x$contact), 1)
  expect_equal(dim(ls_x$sp_contact), c(ls_x$n_gene, ls_x$n_enh), ignore_attr=T)
  expect_equal(as.matrix(ls_x$sp_contact), ls_x$contact, ignore_attr=T)
  expect_equal(dim(ls_x$sp_distance), c(ls_x$n_gene, ls_x$n_enh), ignore_attr=T)
  expect_lt(max(ls_x$sp_distance), 3000000)
})


test_that("data example chr22 no burst sizes", {
  # Options
  opts <- get_opts()
  # No burst sizes
  opts$burst_sizes <- NULL

  chr <- opts$chr_regions
  region <- NULL
  label <- opts$chr_regions
  rlen_unit=1000

  # Tested function: load_data
  ls_x <- IVEA::load_data(opts, chr, region)

  # TEST
  expect_equal(ls_x$n_gene, 193, ignore_attr=T)
  expect_equal(ls_x$n_enh, 3174, ignore_attr=T)
  expect_equal(ls_x$b_size, rep(1, ls_x$n_gene), ignore_attr=T)
})

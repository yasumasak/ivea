
library(optparse)

cat("Run IVEA...\n")

optlist <- list(
  make_option('--chr_regions', action="store", default=NULL, metavar="character", help="Comma-sparated list of chromosomal regions to be analysed. Genes whose TSSs are within this region are analysed. 'chr*:****-****' (start-end) or 'chr*' (whole region)"),
  # Hyperparamters
  make_option('--alpha_pro', action="store", default=80, metavar="numeric", help="Hyperparameter alpha for the gamma prior of promoter activity. [default: %default]"),
  make_option('--alpha_enh', action="store", default=10, metavar="numeric", help="Hyperparameter alpha for the gamma prior of enhancer activity. [default: %default]"),
  make_option('--alpha_o', action="store", default=0.01, metavar="numeric", help="Hyperparameter alpha for the gamma prior of openness. [default: %default]"),
  make_option('--beta_o', action="store", default=0.0001, metavar="numeric", help="Hyperparameter beta for the gamma prior of openness. [default: %default]"),
  make_option('--alpha_k', action="store", default=0.001, metavar="numeric", help="Hyperparameter alpha for the gamma prior of skaler k. [default: %default]"),
  make_option('--beta_k', action="store", default=0.001, metavar="numeric", help="Hyperparameter beta for the gamma prior of skaler k. [default: %default]"),

  # Iteration
  make_option('--max_iteration', action="store", default=500, metavar="numeric", help="Maximum iterations in variational inference. [default: %default]"),
  make_option('--stop_criterion', action="store", default=0.005, metavar="numeric", help="Stop iteration when all relative change of enhancer activites are less than this value. [default: %default]"),

  # Chromatin contact
  make_option('--contacts', action="store", default=NULL, metavar="character", help="File of enhancer-gene contacts. Formatted as 'enhancer-gene_contacts.txt' file produced by 'get_contact_frequency.py'"),
  make_option('--hic_contact_colname', action="store", default=NULL, metavar="character", help="Contact frequency column name in enhancer-gene contacts file. Last column will be analysed if not specified."),
  make_option('--window', action="store", default=5000000, metavar="numeric", help="Analyze enhancer-gene pairs with in this distance. [default: %default]"),

  # Gene expression
  make_option('--rnas', action="store", default=NULL, metavar="character", help="File of Gene expressions."),
  make_option('--rna_id_colname', action="store", default="gene_name", metavar="character", help="Gene id column name in gene expressions file. Gene id must be unique in the file. [default: %default]"),
  make_option('--rna_length_colname', action="store", default="effective_length", metavar="character", help="Effective length column name in gene expressions file. [default: %default]"),
  make_option('--rna_count_colname', action="store", default="expected_count", metavar="character", help="Read count column name in gene expressions file. [default: %default]"),
  make_option('--rna_cutoff_colname', action="store", default="TPM", metavar="character", help="Column name of values examined for expression cutoff in gene expressions file. [default: %default]"),
  make_option('--rna_cutoff_value', action="store", default=8, metavar="numeric", help="Expression cutoff. Genes whose exprssion greater than this value will be analysed. [default: %default]"),

  # Regulatory element
  make_option('--enhancers', action="store", default=NULL, metavar="character", help="File of enhancer elements. Formatted as 'enhancer_elements.txt' produced by 'annotate_regulatory_elements.py'"),
  make_option('--promoters', action="store", default=NULL, metavar="character", help="File of genes promoter elements. Formatted as 'gene_promoter_elements.txt' produced by 'annotate_regulatory_elements.py'"),
  make_option('--burst_sizes', action="store", default=NULL, metavar="character", help="(Optional) Tab-delimited file of gene name (1st column) and Log10(genomic sequence-based burst size) (2nd column). "),
  make_option('--no_le_enh', action="store", default=TRUE, metavar="logical", help="If TRUE, the relative enhancer lengths are NOT used in the prior information for enhancer activities, performing IVEA_nolE. [default: %default]"),

  # Other
  make_option('--outdir', action="store", default=NULL, metavar="character", help="Output directory"),
  make_option('--verbose', action="store_true", default=FALSE, help="If specified, processed data and results in every iterations are saved as RDS in output directory")
)

opts <- parse_args(OptionParser(option_list=optlist), args=commandArgs(trailing=TRUE))

for(p in names(opts)){
  cat("  ",p," = ",opts[[p]],"\n",sep="")
}

# Check required options
if(is.null(opts$chr_regions))  stop("Option 'chr_regions' must be provided.")
if(is.null(opts$contacts))  stop("Option 'contacts' must be provided.")
if(is.null(opts$rnas))      stop("Option 'rnas' must be provided.")
if(is.null(opts$enhancers)) stop("Option 'enhancers' must be provided.")
if(is.null(opts$promoters)) stop("Option 'promoters' must be provided.")
if(is.null(opts$outdir)) stop("Option 'outdir' must be provided.")

# Output directory
dir.create(opts$outdir, showWarnings=F, recursive=T)

# Parse chromosomal region to be analysed
parse_chr_region <- function(chr_region){
  v_cr <- unlist(strsplit(chr_region,":"))
  if(length(v_cr)==1){
    list(chr=v_cr[1], region=NULL, label=v_cr[1])
  }else{
    list(chr=v_cr[1], region=as.numeric(unlist(strsplit(v_cr[2],"-"))), label=paste0(v_cr[1],".",v_cr[2]))
  }
}


#------------------------------
# Main
#------------------------------
v_chr_region <- unlist(strsplit(opts$chr_regions,","))
for(chr_region in v_chr_region){
  cat("\n  For",chr_region,"\n")
  t <- proc.time()

  # Chromosomal region to be analysed
  cr <- parse_chr_region(chr_region)
  chr <- cr$chr
  region <- cr$region
  label <- cr$label

  #------------------------------
  # Load data
  ls_x <- IVEA::load_data(opts, chr, region)
  if(opts$verbose){
    # Save the data
    rds_data <- paste0(opts$outdir,"/ivea_data.",label,".rds")
    saveRDS(ls_x, file=rds_data)
  }

  #------------------------------
  # Run variational Bayes
  ls_z <- IVEA::run_variational_bayes(opts, ls_x)
  if(opts$verbose){
    # Save the estimates
    rds_estimates <- paste0(opts$outdir,"/ivea_estimates.",label,".rds")
    saveRDS(ls_z, file=rds_estimates)
  }

  #------------------------------
  # Output
  ls_out <- IVEA::make_summary(ls_x, ls_z)
  out_info <- paste0(opts$outdir,"/predictions_info.",label,".txt")
  out_bedpe <- paste0(opts$outdir,"/predictions_score.",label,".bedpe")
  out_bed_pro_bed <- paste0(opts$outdir,"/estimates.promoter_activity.",label,".bed")
  out_bed_enh_bed <- paste0(opts$outdir,"/estimates.enhancer_activity.",label,".bed")
  out_bed_pro_all <- paste0(opts$outdir,"/estimates.promoter_activity.",label,".txt")
  out_bed_enh_all <- paste0(opts$outdir,"/estimates.enhancer_activity.",label,".txt")
  write.table(ls_out$info_pair, file=out_info, sep="\t", col.names=T, row.names=F, quote=F)
  write.table(ls_out$bedpe_pair, file=out_bedpe, sep="\t", col.names=F, row.names=F, quote=F)
  write.table(ls_out$bed_pro[, 1:5], file=out_bed_pro_bed, sep="\t", col.names=F, row.names=F, quote=F)
  write.table(ls_out$bed_enh[, 1:5], file=out_bed_enh_bed, sep="\t", col.names=F, row.names=F, quote=F)
  write.table(ls_out$bed_pro, file=out_bed_pro_all, sep="\t", col.names=T, row.names=F, quote=F)
  write.table(ls_out$bed_enh, file=out_bed_enh_all, sep="\t", col.names=T, row.names=F, quote=F)

  cat("  Elapsed time: ",(proc.time() - t)[3],"\n", sep="")
}
cat("\nDone.\n\n")


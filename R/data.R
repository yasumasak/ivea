
#' Load data files
#'
#' Load data files according to configuration
#'
#' @param config a list object of configuration containing the following components:
#'   \describe{
#'     \item{contacts}{file of Enhancer-gene contacts.}
#'     \item{hic_contact_colname}{contact frequency column name in enhancer-gene
#'       contacts file.}
#'     \item{window}{analyze enhancer-gene pairs with in this distance.}
#'     \item{rnas}{file of gene expressions.}
#'     \item{rna_id_colname}{gene id column name in gene expressions file.
#'       Gene id must be unique in the file.}
#'     \item{rna_length_colname}{transcript length column name in gene expressions
#'       file.}
#'     \item{rna_count_colname}{read count column name in gene expressions file.}
#'     \item{rna_cutoff_colname}{column name of values examined for expression
#'       cutoff in gene expressions file.}
#'     \item{rna_cutoff_value}{expression cutoff. Genes whose exprssion greater
#'       than this value will be analysed.}
#'     \item{enhancers}{file of enhancer elements.}
#'     \item{promoters}{file of genes promoter elements.}
#'     \item{burst_sizes}{file of log10(genomic sequence-based burst size).}
#'   }
#' @param chr chromosome to be analyzed
#' @param region region on chromosome to be analyzed
#' @param rlen_unit unit length used for transcript, enhancer and promoter length
#'
#' @return A list object of input data for the variational inference,
#'   containing the following components:
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
#' @export
load_data <- function(config, chr, region, rlen_unit=1000){
  cat("  Loading data\n")

  # Genomic sequence-based burst sizes
  if(!is.null(config$burst_sizes)){
    df_log10_bs <- data.table::fread(config$burst_sizes, header=F, sep="\t",
                                     stringsAsFactors=F, check.names=F, quote="", data.table=F)
    v_log10_bs <- df_log10_bs[, 2]
    names(v_log10_bs) <- df_log10_bs[, 1]
  }else{
    v_log10_bs <- NULL
  }

  # RNA-seq
  df_rna <- data.table::fread(config$rnas, header=T, sep="\t",
                              stringsAsFactors=F, check.names=F, quote="", data.table=F)
  # Gene promoter elements
  df_pro <- data.table::fread(config$promoters, header=T, sep="\t",
                              stringsAsFactors=F, check.names=F, quote="", data.table=F)
  # Gene body lengths
  df_enh <- data.table::fread(config$enhancers, header=T, sep="\t",
                              stringsAsFactors=F, check.names=F, quote="", data.table=F)
  # Contact frequency
  df_con_all <- data.table::fread(config$contacts, header=T, sep="\t",
                                  stringsAsFactors=F, check.names=F, quote="", data.table=F)

  rownames(df_rna) <- df_rna[, config$rna_id_colname]
  rownames(df_pro) <- df_pro$gene_name
  rownames(df_enh) <- df_enh$name

  # Gene TSS
  df_gene_tss <- df_pro[, c('gene_name', 'tss')]

  # Set contact column
  if(is.null(config$hic_contact_colname)){
    contact_colname <- colnames(df_con_all)[ncol(df_con_all)] # Using the last column
  }else{
    contact_colname <- config$hic_contact_colname # Using specified contact column
  }
  cat("  Contact frequency column name:",contact_colname,"\n")
  df_con_all <- df_con_all[, c('chr', 'start', 'end', 'name', 'target_gene',
                                 'target_tss', 'distance', contact_colname)]
  colnames(df_con_all)[8] <- 'contact'

  #------------------------------
  # Genes to be analysed
  #------------------------------
  # Extract contacts by specified chromosomal region and distance cutoff
  if(is.null(region)){
    df_con_lim <- df_con_all[df_con_all$chr == chr & df_con_all$distance < config$window, ]
  }else{
    df_con_lim <- df_con_all[df_con_all$chr == chr & df_con_all$distance < config$window
                                                  & df_con_all$target_tss > region[1]
                                                  & df_con_all$target_tss < region[2], ]
  }

  # Gene list
  v_g <- unique(df_con_lim[, 'target_gene'])
  cat("  Genes in the region of interest:",length(v_g),"\n")

  # Genes expressed
  if(is.na(config$rna_cutoff_value)){
    df_expressed <- df_rna
  }else{
    df_expressed <- df_rna[df_rna[, config$rna_cutoff_colname] > config$rna_cutoff_value, ]
  }
  v_g_expressed <- df_expressed[stats::na.omit(match(v_g, df_expressed[, config$rna_id_colname])),
                                config$rna_id_colname] # omit genes not expressed
  v_g_omit <- v_g[!v_g %in% v_g_expressed]
  cat("  Genes with insufficient expression:",length(v_g_omit),"\n")

  # Genes to be analyzed
  v_gene <- v_g_expressed

  #------------------------------
  # Re-extract data based on
  # the genes to be analysed
  #------------------------------
  # Gene TSS
  v_tss <- df_gene_tss[v_gene, 'tss']

  # Burst sizes
  if(!is.null(v_log10_bs)){
    v_b_size <- 10 ^ v_log10_bs[v_gene]
  }else{
    v_b_size <- rep(1, length(v_gene))
    names(v_b_size) <- v_gene
  }

  # RNA
  v_rna_count <- round(df_rna[v_gene, config$rna_count_colname])
  v_rna_length <- df_rna[v_gene, config$rna_length_colname]

  # Promoter elements
  v_pro <- df_pro[v_gene, 'name']
  v_pro_count <- df_pro[v_gene, 'read_count']
  v_pro_length <- df_pro[v_gene, 'length']

  # Enhancer elements
  df_con <- df_con_lim[df_con_lim$target_gene %in% v_gene, ]
  v_enh <- unique(df_con[, 'name'])
  v_enh_count <- df_enh[v_enh, 'read_count']
  v_enh_length <- df_enh[v_enh, 'length']

  # Contact frequency
  f.gene <- factor(df_con$target_gene, levels=v_gene)
  f.enh <- factor(df_con$name, levels=v_enh)
  v_i <- as.numeric(f.gene)
  v_j <- as.numeric(f.enh)
  sp_dist <- Matrix::sparseMatrix(i=v_i, j=v_j, x=df_con$distance)
  sp_con <- Matrix::sparseMatrix(i=v_i, j=v_j, x=df_con$contact)
  mx_con <- as.matrix(sp_con)

  # Relative length for effective rna length
  v_rna_rlen <- v_rna_length / rlen_unit
  # Relative length for regulatory elements
  v_pro_rlen <- v_pro_length / rlen_unit
  v_enh_rlen <- v_enh_length / rlen_unit

  # Number of genes and regulatory elements
  n_gene <- length(v_gene)
  n_enh <- length(v_enh)
  n_pro <- length(v_pro)

  if(n_gene != n_pro){stop("The numbers of genes and promters are different.")}
  cat("  Genes to be analysed:",n_gene,"\n")
  cat("  Enhancers to be analysed:",n_enh,"\n")

  # Input data for IVEA
  ls_x <- list(chr=chr, region=region,
               gene_name=v_gene, gene_tss=v_tss,
               rna_count=v_rna_count, rna_rlen=v_rna_rlen,
               pro_name=v_pro, pro_count=v_pro_count, pro_rlen=v_pro_rlen,
               enh_name=v_enh, enh_count=v_enh_count, enh_rlen=v_enh_rlen,
               n_gene=n_gene, n_enh=n_enh,
               b_size=v_b_size,
               contact=mx_con,
               sp_contact=sp_con, sp_distance=sp_dist)

  return(ls_x)
}

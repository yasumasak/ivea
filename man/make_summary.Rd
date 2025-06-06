% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IVEA.R
\name{make_summary}
\alias{make_summary}
\title{Make summary data frame}
\usage{
make_summary(ls_x, ls_z)
}
\arguments{
\item{ls_x}{an input-data list object containing the following components:
\describe{
\item{chr}{chromosome.}
\item{region}{region on chromosome.}
\item{gene_name}{character vector of gene names.}
\item{gene_tss}{numeric vector of gene TSSs.}
\item{rna_count}{numeric vector of sequence read counts on genes.}
\item{rna_rlen}{numeric vector of relative effective lengths of genes.}
\item{pro_name}{character vector of gene promoter names.}
\item{pro_count}{numeric vecor of sequence read counts gene promoters.}
\item{pro_rlen}{numeric vecor of relative lengths of gene promoters.}
\item{enh_name}{character vector of enhancer names.}
\item{enh_count}{numeric vecor of sequence read counts enhancers.}
\item{enh_rlen}{numeric vecor of relative lengths of enhancers.}
\item{n_gene}{total number of genes.}
\item{n_enh}{total number of enhancers.}
\item{b_size}{numeric vector of genomic sequence-based burst size.}
\item{contact}{matrix (genes x enhancers) of contact frequencies.}
\item{sp_contact}{sparse matrix (genes x enhancers) of contact frequencies.}
\item{sp_distance}{sparse matrix (genes x enhancers) of genomic distances.}
}}

\item{ls_z}{a resultant list object containing the following components:
\describe{
\item{k}{vector of estimated scaler variable.}
\item{pro_openness}{matrix of estimated opennesses of gene promoters.}
\item{enh_openness}{matrix of estimated opennesses of enhancers.}
\item{pro_inv_openness}{matrix of estimated 1/opennesses of gene promoters.}
\item{enh_inv_openness}{matrix of estimated 1/opennesses of enhancers.}
\item{pro_activity}{matrix of estimated activities of gene promoters.}
\item{enh_activity}{matrix of estimated activities of enhancers.}
\item{pro_log_activity}{matrix of estimated log activities of gene promoters.}
\item{enh_log_activity}{matrix of estimated log activities of enhancers.}
}}
}
\value{
A list object containing the following summary data frames:
\describe{
\item{info_pair}{data frame of information of enhancer-gene pairs.}
\item{bedpe_pair}{data frame of prediction scores in BEDPE format.}
}
}
\description{
Calculate enhancer-gene pairs' strength, contribution and score,
and make summary data frames of enhancer-gene pairs
}

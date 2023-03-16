import argparse
import pandas as pd
import numpy as np
import os
import os.path
from tools import *


def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description="Get regulatory elements.", epilog=epilog, formatter_class=formatter)
    readable = argparse.FileType('r')

    parser.add_argument('--outdir', required=required_args, help="Directory to write output files to.")

    # Inputs
    parser.add_argument('--narrowPeak', required=required_args, help="narrowPeak file output by macs2. Must include summits (--call-summits)")
    parser.add_argument('--bam', required=required_args, help="Bam file of DNAase-Seq or ATAC-Seq")
    parser.add_argument('--chrom_sizes', required=required_args, help="File listing chromosome sizes")
    parser.add_argument('--regions_whitelist', default="", help="Bed file of regions to forcibly include in regulatory elements. Overrides regions_blacklist")
    parser.add_argument('--regions_blacklist', default="", help="Bed file of regions to forcibly exclude from regulatory elements")

    # Parameters for peaks
    parser.add_argument('--n_peaks', default=150000, help="Number of peaks to use for defining regulatory elements")
    parser.add_argument('--extension_from_summit', default=250, help="Number of base pairs to extend each preak from its summit (or from both ends of region if using --ignore_summit)")
    parser.add_argument('--ignore_summit', action="store_true", help="Compute peaks using the full peak regions, rather than extending from summit.")
    parser.add_argument('--min_peak_width', default=500, help="Candidate elements whose width is below this threshold are expanded to this width. Only used with --ignore_summit")

    # Gene annotations
    parser.add_argument('--genes', required=required_args, help="Bed-6 file with gene annotations. Will be used to assign TSS to genes.")
    parser.add_argument('--gene_annotations', default="symbol", help="Comma delimited string of names corresponding to what present in the 'name' field of the gene annotation bed file")
    parser.add_argument('--gene_identifier', default="symbol", help="Identifier used to identify genes. Must be present in gene_annotations. The gene identifier must be unique")

    # Other
    parser.add_argument('--tss_slop_for_class_assignment', default=1000, type=int, help="Consider an element a promoter if it is within this many bp of a tss")
    parser.add_argument('--use_secondary_counting_method', action="store_true", help="Use a slightly slower way to count bam over bed. Also requires more memory. But is more stable")

    args = parser.parse_args()
    return(args)


def main(args):

    print("Get regulatory elements...")

    os.makedirs(args.outdir, exist_ok=True)
    print_params(args)
    print("\n")

    # Make regulatory elements
    if not args.ignore_summit:
        elements_bed = make_regulatory_elements_from_summits(macs_peaks = args.narrowPeak,
                                            bam_file = args.bam,
                                            chrom_sizes = args.chrom_sizes,
                                            regions_whitelist = args.regions_whitelist,
                                            regions_blacklist = args.regions_blacklist,
                                            n_peaks = args.n_peaks,
                                            peak_extend = args.extension_from_summit,
                                            outdir = args.outdir)
    else:
        elements_bed = make_regulatory_elements_from_peaks(macs_peaks = args.narrowPeak,
                                    bam_file = args.bam,
                                    chrom_sizes = args.chrom_sizes,
                                    regions_whitelist = args.regions_whitelist,
                                    regions_blacklist = args.regions_blacklist,
                                    n_peaks = args.n_peaks,
                                    peak_extend = args.extension_from_summit,
                                    min_peak_width = args.min_peak_width,
                                    outdir = args.outdir)

    # Setup genes
    genes = load_genes(genes_bed = args.genes,
                       chrom_sizes = args.chrom_sizes,
                       outdir = args.outdir,
                       gene_annotations = args.gene_annotations,
                       gene_identifier = args.gene_identifier)

    # Setup Regulatory elements
    annotate_elements(genes=genes,
                  chrom_sizes=args.chrom_sizes,
                  bam_file=args.bam,
                  elements_bed=elements_bed,
                  tss_slop_for_class_assignment=args.tss_slop_for_class_assignment,
                  use_fast_count = (not args.use_secondary_counting_method),
                  outdir = args.outdir)

    print("Done.\n")


def make_regulatory_elements_from_summits(macs_peaks, bam_file, chrom_sizes, regions_whitelist, regions_blacklist, n_peaks, peak_extend, outdir):
    ## Generate regulatory elements from MACS summits

    sorted_macs_peaks = os.path.join(outdir, os.path.basename(macs_peaks) + ".sorted")
    counts_bedgraph = os.path.join(outdir, os.path.basename(macs_peaks) + ".read_counts.bedgraph")
    elements_bed = os.path.join(outdir, "regulatory_elements.bed")

    if regions_whitelist:
        whitelist_command = "(bedtools intersect -a {regions_whitelist} -b {chrom_sizes}.bed -wa | cut -f 1-3 && cat) |"
    else:
        whitelist_command = ""

    if regions_blacklist:
        blacklist_command = "bedtools intersect -v -wa -a stdin -b {regions_blacklist} | "
    else:
        blacklist_command = ""

    # 0. Sort DHS/ATAC reads in peaks
    sort_bed(macs_peaks, sorted_macs_peaks, chrom_sizes)

    # 1. Count DHS/ATAC reads in peaks
    run_count_reads(bam_file, counts_bedgraph, sorted_macs_peaks, chrom_sizes, use_fast_count=True)

    # 2. Take top N peaks, get summits, extend summits, merge, remove blacklist, add whitelist, sort and merge
    # use -sorted in intersect command? Not worth it, both files are small
    command = "bedtools sort -i {counts_bedgraph} -faidx {chrom_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_peaks} | " + \
        "bedtools intersect -b stdin -a {sorted_macs_peaks} -wa | " + \
        "awk '{{print $1 \"\\t\" $2 + $10 \"\\t\" $2 + $10}}' | " + \
        "bedtools slop -i stdin -b {peak_extend} -g {chrom_sizes} | " + \
        "bedtools sort -i stdin -faidx {chrom_sizes} | " + \
        "bedtools merge -i stdin | " + \
        blacklist_command + \
        "cut -f 1-3 | " + whitelist_command + \
        "bedtools sort -i stdin -faidx {chrom_sizes} | bedtools merge -i stdin > {elements_bed}"

    command = command.format(**locals())

    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("  Running: " + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')
    if not err == '':
        raise RuntimeError("Command failed.")

    return elements_bed


def make_regulatory_elements_from_peaks(macs_peaks, bam_file, chrom_sizes, regions_whitelist, regions_blacklist, n_peaks, peak_extend, min_peak_width, outdir):
    ## Generate regulatory elements from MACS narrowPeak - do not use summits

    sorted_macs_peaks = os.path.join(outdir, os.path.basename(macs_peaks) + ".sorted")
    counts_bedgraph = os.path.join(outdir, os.path.basename(macs_peaks) + ".read_counts.bedgraph")
    elements_bed = os.path.join(outdir, "regulatory_elements.bed")

    if regions_whitelist:
        whitelist_command = "(bedtools intersect -a {regions_whitelist} -b {chrom_sizes}.bed -wa | cut -f 1-3 && cat) |"
    else:
        whitelist_command = ""

    if regions_blacklist:
        blacklist_command = "bedtools intersect -v -wa -a stdin -b {regions_blacklist} | "
    else:
        blacklist_command = ""

    # 0. Sort DHS/ATAC reads in peaks
    sort_bed(macs_peaks, sorted_macs_peaks, chrom_sizes)

    # 1. Count DHS/ATAC reads in peaks
    run_count_reads(bam_file, counts_bedgraph, sorted_macs_peaks, chrom_sizes, use_fast_count=True)

    # 2. Take top N peaks, extend peaks (min size 500), merge, remove blacklist, add whitelist, sort and merge
    # use -sorted in intersect command? Not worth it, both files are small
    command = "bedtools sort -i {counts_bedgraph} -faidx {chrom_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_peaks} | " + \
        "bedtools intersect -b stdin -a {sorted_macs_peaks} -wa | " + \
        "bedtools slop -i stdin -b {peak_extend} -g {chrom_sizes} | " + \
        "awk '{{ l=$3-$2; if (l < {min_peak_width}) {{ $2 = $2 - int(({min_peak_width}-l)/2); $3 = $3 + int(({min_peak_width}-l)/2) }} print $1 \"\\t\" $2 \"\\t\" $3}}' | " + \
        "bedtools sort -i stdin -faidx {chrom_sizes} | " + \
        "bedtools merge -i stdin | " + \
        blacklist_command + \
        "cut -f 1-3 | " + whitelist_command + \
        "bedtools sort -i stdin -faidx {chrom_sizes} | bedtools merge -i stdin > {elements_bed}"

    command = command.format(**locals())

    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("  Running: " + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')
    if not err == '':
        raise RuntimeError("Command failed.")

    return elements_bed


def load_genes(genes_bed,
               chrom_sizes,
               outdir,
               gene_annotations,
               gene_identifier):

    bed = read_bed6(genes_bed)
    genes = process_gene_bed(bed, gene_annotations, gene_identifier, chrom_sizes)

    return genes


def process_gene_bed(bed, annotations, identifier, chrom_sizes=None, fail_on_nonunique=True):

    assert(identifier in annotations)

    names = bed.name.str.split(";", expand=True)
    assert(len(names.columns) == len(annotations.split(",")))
    names.columns = annotations.split(",")
    bed = pd.concat([bed, names], axis=1)

    bed.rename(columns={'name': 'gene_name'}, inplace=True)
    bed['gene_name'] = bed[identifier]
    bed['tss'] = get_tss_for_bed(bed)

    bed.drop_duplicates(inplace=True)

    # Remove genes that are not defined in chromosomes file
    if chrom_sizes is not None:
        sizes = read_bed6(chrom_sizes)
        bed['chr'] = bed['chr'].astype('str')
        bed = bed[bed['chr'].isin(set(sizes['chr'].values))]

    # Enforce that gene names should be unique
    if fail_on_nonunique:
        assert(len(set(bed['gene_name'])) == len(bed['gene_name'])), \
               "Gene IDs are not unique! Failing. Please ensure unique identifiers are passed to --genes"

    return bed


def annotate_elements(outdir=".",
                  chrom_sizes="",
                  bam_file="",
                  genes=None,
                  elements_bed="",
                  tss_slop_for_class_assignment = 1000,
                  force=True,
                  use_fast_count=True):

    elements = read_bed6(elements_bed)
    elements['chr'] = elements['chr'].astype('str')

    elements['length'] = abs(elements['end'] - elements['start'])
    elements = add_counts_for_bed(elements, elements_bed, chrom_sizes, bam_file, outdir, "regulatory_elements", force, use_fast_count)

    # Assign categories
    if genes is not None:
        print("  Assigning classes to elements")
        elements = assign_element_classes(elements, genes, tss_slop = tss_slop_for_class_assignment)
        promoters = elements.loc[elements['class'] == "promoter", :].copy()

    # Output enhancer elements
    elements.to_csv(os.path.join(outdir, "enhancer_elements.txt"), sep="\t", index=False, header=True, float_format="%.6f")
    elements[['chr','start','end','name']].to_csv(os.path.join(outdir, "enhancer_elements.bed"), sep="\t", index=False, header=False)

    # Output promoter elements
    if genes is not None:
        gene_promoters = get_gene_promoter(promoters, genes, tss_slop = tss_slop_for_class_assignment)
        gene_promoters.to_csv(os.path.join(outdir, "gene_promoter_elements.txt"), sep="\t", index=False, header=True, float_format="%.6f")


def assign_element_classes(elements, genes, tss_slop=1000):

    # build pyranges df
    pr_tss = df_to_pyranges(genes, start_col='tss', end_col='tss', start_slop=tss_slop, end_slop=tss_slop)
    pr_gene = df_to_pyranges(genes)

    # label everything as intergenic
    elements['class'] = "intergenic"
    elements['uid'] = range(elements.shape[0])
    pr_enh = df_to_pyranges(elements)

    # genic element
    pr_genic_enh = pr_enh.join(pr_gene, suffix="_genic")
    df_genic_enh = pr_genic_enh.df[['gene_name','uid']].groupby('uid',as_index=False).aggregate(lambda x: ','.join(list(set(x))))

    # promoter element
    pr_promoter_enh = pr_enh.join(pr_tss, suffix="_promoter")
    df_promoter_enh = pr_promoter_enh.df[['gene_name','uid']].groupby('uid',as_index=False).aggregate(lambda x: ','.join(list(set(x))))

    # set class
    elements.loc[elements['uid'].isin(df_genic_enh.uid), 'class'] = "genic"
    elements.loc[elements['uid'].isin(df_promoter_enh.uid), 'class'] = "promoter"

    # output stats
    print("  Total elements: {}".format(len(elements)))
    print("        Promoters: {}".format(sum(elements['class'] == "promoter")))
    print("        Genic: {}".format(sum(elements['class'] == "genic")))
    print("        Intergenic: {}".format(sum(elements['class'] == "intergenic")))

    # add promoter/genic gene
    elements = elements.merge(df_promoter_enh.rename(columns={'gene_name':'promoter_gene'}), on='uid', how='left').fillna(value={'promoter_gene':""})
    elements = elements.merge(df_genic_enh.rename(columns={'gene_name':'genic_gene'}), on='uid', how='left').fillna(value={'genic_gene':""})
    elements.drop(['uid'], axis=1, inplace=True)

    # just to keep things consistent with original code 
    elements['name'] = elements.apply(lambda e: "{}:{}-{}".format(e.chr, e.start, e.end), axis=1)
    return elements


def get_gene_promoter(promoters, genes, tss_slop=1000):

    # build pyranges df
    pr_tss = df_to_pyranges(genes, start_col='tss', end_col='tss', start_slop=tss_slop, end_slop=tss_slop)
    pr_pro = df_to_pyranges(promoters)

    # promoters in each gene
    pr_gene_pro = pr_pro.join(pr_tss, suffix='_gene')
    df_gene_pro = pr_gene_pro.df[['chr','tss','strand','gene_name','name','length','read_count']].groupby('gene_name', as_index=False)\
                  .aggregate({'chr':['first'], 'tss':['first'], 'strand':['first'], 'name':lambda x: ','.join(list(set(x))), 'length':['sum'], 'read_count':['sum']})
    gene_promoters = df_gene_pro.set_axis(['gene_name','chr','tss','strand','name','length','read_count'], axis='columns')

    return gene_promoters


if __name__ == '__main__':
    args = parseargs()
    main(args)


import argparse
import pandas as pd
import os, os.path
from gtfparse import read_gtf


def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description="Make gencode annotation bed file.", epilog=epilog, formatter_class=formatter)
    readable = argparse.FileType('r')

    parser.add_argument('--gtf_gencode', required=required_args, help="GTF annotation file of Gencode.")
    parser.add_argument('--outdir', required=True, help="Directory to write output files to.")
    parser.add_argument('--chromosomes', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX", help="Comma-sparated list of chromosomes included in the output.")

    args = parser.parse_args()
    return(args)


def main(args):

    os.makedirs(args.outdir, exist_ok=True)

    # Load Gencode annotation
    df_gencode = read_gtf(args.gtf_gencode)
    # Types to use
    in_type = ['bidirectional_promoter_lncRNA', '3prime_overlapping_ncRNA', 'polymorphic_pseudogene', 'transcribed_unitary_pseudogene', 'TEC', \
        'unitary_pseudogene', 'sense_overlapping', 'transcribed_processed_pseudogene', 'processed_transcript', 'pseudogene', 'sense_intronic', \
        'transcribed_unprocessed_pseudogene', 'unprocessed_pseudogene', 'antisense', 'lincRNA', 'processed_pseudogene', 'protein_coding']

    # Split chromosome list
    chromosomes = args.chromosomes.split(",")
    list_all_gene_id = []
    list_all_gene_name = []
    list_all_id2name = []

    for chromosome in chromosomes:
        # gene_id based
        df_gene_id = df_gencode.loc[(df_gencode['feature'] == "gene") & (df_gencode['gene_type'].isin(in_type)) & (df_gencode['seqname'] == chromosome), \
                                   ['seqname', 'start', 'end', 'gene_id', 'gene_type', 'strand']]
        df_gene_id['start'] = df_gene_id['start'] - 1 # Change 1-based to 0-based positions
        list_all_gene_id.append(df_gene_id)

        # gene_name based
        df_gene_name = df_gencode.loc[(df_gencode['feature'] == "gene") & (df_gencode['gene_type'].isin(in_type)) & (df_gencode['seqname'] == chromosome), \
                                     ['seqname', 'start', 'end', 'gene_name', 'gene_type', 'strand']]
        df_gene_name['start'] = df_gene_name['start'] - 1 # Change 1-based to 0-based positions
        list_all_gene_name.append(df_gene_name)

        # list of gene_id to gene_name
        df_id2name = df_gencode.loc[(df_gencode['feature'] == "gene") & (df_gencode['gene_type'].isin(in_type)) & (df_gencode['seqname'] == chromosome), \
                                     ['gene_id', 'gene_name']]
        list_all_id2name.append(df_id2name)

    df_all_gene_id = pd.concat(list_all_gene_id)
    df_all_gene_name = pd.concat(list_all_gene_name)
    df_all_id2name = pd.concat(list_all_id2name)

    # Output file
    df_all_gene_id.to_csv(os.path.join(args.outdir, "gencode.annotation.gene_id.bed"), sep="\t", index=False, header=False, float_format="%.6f")
    df_all_gene_name.drop_duplicates(subset='gene_name', keep=False)\
      .to_csv(os.path.join(args.outdir, "gencode.annotation.gene_name.bed"), sep="\t", index=False, header=False, float_format="%.6f")
    df_all_id2name.to_csv(os.path.join(args.outdir, "gencode.gene_id2name.txt"), sep="\t", index=False, header=False, float_format="%.6f")


if __name__ == '__main__':
    args = parseargs()
    main(args)


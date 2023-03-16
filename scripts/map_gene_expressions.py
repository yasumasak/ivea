import argparse
import pandas as pd
import numpy as np
import os, os.path
from gtfparse import read_gtf
from tools import *

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description="Make gene expression information file.", epilog=epilog, formatter_class=formatter)
    readable = argparse.FileType('r')

    parser.add_argument('--bed_ref', required=required_args, help="Bed-6 file of referential gene annotations.")
    parser.add_argument('--gtf_gencode', required=required_args, help="GTF from Gencode. Should be the same one used for RSEM.")
    parser.add_argument('--rsem', required=required_args, help="genes.results from RSEM.")
    parser.add_argument('--outdir', required=True, help="Directory to write output files to.")

    args = parser.parse_args()
    return(args)


def load_gencode_rsem(gtf_gencode, rsem):

    # RSEM genes.results
    df_rsem = pd.read_csv(rsem, sep='\t', header=0).loc[:, ['gene_id', 'effective_length', 'expected_count', 'TPM']]

    # Gencode GTF file
    df_gencode = read_gtf(gtf_gencode)
    in_type = ['bidirectional_promoter_lncRNA', '3prime_overlapping_ncRNA', 'polymorphic_pseudogene', 'transcribed_unitary_pseudogene', 'TEC', \
        'unitary_pseudogene', 'sense_overlapping', 'transcribed_processed_pseudogene', 'processed_transcript', 'pseudogene', 'sense_intronic', \
        'transcribed_unprocessed_pseudogene', 'unprocessed_pseudogene', 'antisense', 'lincRNA', 'processed_pseudogene', 'protein_coding']
    df_gene = df_gencode.loc[(df_gencode['feature'] == 'gene') & (df_gencode['gene_type'].isin(in_type)) & (df_gencode['seqname'] != 'chrY'), \
                             ['seqname', 'start', 'end', 'gene_id', 'strand', 'gene_name', 'gene_type']]
    df_gene['start'] = df_gene['start'] - 1 # Change 1-based to 0-based positions

    # Genes with RSEM information
    df_gene_rsem = pd.merge(df_gene, df_rsem, how='inner', on='gene_id')

    return(df_gene_rsem)


def main(args):

    print("Get gene expressions...")

    os.makedirs(args.outdir, exist_ok=True)
    print_params(args)
    print("\n")

    # Load Gencode information and RSEM result
    df_gene_rsem = load_gencode_rsem(args.gtf_gencode, args.rsem)

    # Load Bed file of referential gene annotations
    df_ref = pd.read_csv(args.bed_ref, sep="\t", header=None)
    df_ref.columns = ['chr', 'start', 'end', 'gene_name', 'score', 'strand']

    # 1. Combine by gene_name --------------------
    df_merged_1 = pd.merge(df_ref, df_gene_rsem.loc[:, ['gene_name', 'gene_id', 'effective_length', 'expected_count', 'TPM']],\
                           how='inner', on='gene_name').drop('score', axis=1)

    # Remained data
    df_ref_r1 = df_ref[np.logical_not(df_ref['gene_name'].isin(df_merged_1['gene_name']))]
    df_gene_rsem_r1 = df_gene_rsem[np.logical_not(df_gene_rsem['gene_name'].isin(df_merged_1['gene_name']))]
    # Pyranges for examining overlaps
    gr_ref_r1 = df_to_pyranges(df_ref_r1, strand_col='strand')
    gr_gene_rsem_r1 = df_to_pyranges(df_gene_rsem_r1, chr_col='seqname', strand_col='strand')
    # Join overlapped records
    df_joined = gr_ref_r1.join(gr_gene_rsem_r1, strandedness='same', report_overlap=True).df
    # Overlapped fractions
    frac_overlap = df_joined['Overlap'] / (df_joined['end'] - df_joined['start'])

    # 2. Combined by overlaps --------------------
    df_merged_2 = df_joined.loc[frac_overlap > 0.9, ['chr', 'start', 'end', 'gene_name', 'strand', 'gene_id', 'effective_length', 'expected_count', 'TPM']]

    # Output file
    pd.concat([df_merged_1, df_merged_2]).drop_duplicates(subset='gene_name', keep=False)\
      .to_csv(os.path.join(args.outdir, "gene_expressions.txt"), sep="\t", index=False, header=True, float_format="%.6f")

if __name__ == '__main__':
    args = parseargs()
    main(args)


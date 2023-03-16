import argparse
import pandas as pd
import numpy as np
import os, os.path
from tools import *

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description="Make gene expression information file.", epilog=epilog, formatter_class=formatter)
    readable = argparse.FileType('r')

    parser.add_argument('--table', required=required_args, help="Table of gene id and gene name.")
    parser.add_argument('--rsem', required=required_args, help="genes.results from RSEM.")
    parser.add_argument('--outdir', required=True, help="Directory to write output files to.")

    args = parser.parse_args()
    return(args)


def main(args):

    print("Get gene expressions...")

    os.makedirs(args.outdir, exist_ok=True)
    print_params(args)
    print("\n")

    # RSEM genes.results
    df_rsem = pd.read_csv(args.rsem, sep='\t', header=0).loc[:, ['gene_id', 'effective_length', 'expected_count', 'TPM']]

    # Load table of gene id - gene name
    df_id2name = pd.read_csv(args.table, sep="\t", header=None)
    df_id2name.columns = ['gene_id', 'gene_name']

    # 1. Combine by gene_name --------------------
    df_merged = pd.merge(df_id2name, df_rsem, how='inner', on='gene_id')

    # Output file
    df_merged.drop_duplicates(subset='gene_name', keep=False)\
      .to_csv(os.path.join(args.outdir, "gene_expressions.txt"), sep="\t", index=False, header=True, float_format="%.6f")

if __name__ == '__main__':
    args = parseargs()
    main(args)


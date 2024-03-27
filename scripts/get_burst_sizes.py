import argparse
import pandas as pd
import numpy as np
import os, os.path
from tools import *

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description="Calculate burst size estimates.", epilog=epilog, formatter_class=formatter)
    readable = argparse.FileType('r')

    parser.add_argument('--epd_bed_file', '-b', required=True, help="Bed file from EPD.")
    parser.add_argument('--epd_motif_file', '-m', required=True, help="promoter_motifs.txt file from EPD.")
    parser.add_argument('--epd_bed_file_2', '-b2', required=False, help="(Optional) Second bed file from EPD.")
    parser.add_argument('--epd_motif_file_2', '-m2', required=False, help="(Optional) Second promoter_motifs.txt file from EPD.")
    parser.add_argument('--genes', required=True, help="Bed-6 file of gene annotations. Should be the same one used in used in get_regulatory_elements.py.")
    parser.add_argument('--outdir', required=True, help="Directory to write output files to.")
    parser.add_argument('--tss_slop', default=1000, type=int, help="Consider a region within this many bp from a tss as a promoter region.")

    args = parser.parse_args()
    return(args)


def load_epd(epd_bed_file, motif_file):

    # EPD promoters
    df_epd = pd.read_csv(epd_bed_file, sep="\t", header=None).iloc[:, 0:6]
    df_epd.columns = ['chr', 'start', 'end', 'epd_id', 'score', 'strand']
    # Motifs in EPD promoters
    df_motif = pd.read_csv(motif_file, sep="\t", header=0)
    df_motif.columns = ['epd_id', 'TATA-box', 'Inr', 'CCAAT-box', 'GC-box']
    # Merge
    epd_motif = pd.merge(df_epd, df_motif, how='left', on='epd_id')

    return(epd_motif)


def get_gene_info(genes):

    # Gene model used in IVEA
    df_gene = read_bed6(genes)
    df_gene.columns = ['chr', 'start', 'end', 'gene_name', 'score', 'strand']
    df_gene['tss'] = get_tss_for_bed(df_gene)
    df_gene['gene_body_length'] = df_gene['end'] - df_gene['start']

    return(df_gene)


def main(args):

    print("Get burst sizes...")

    os.makedirs(args.outdir, exist_ok=True)
    print_params(args)
    print("\n")

    # Load EPD files
    epd_motif_1 = load_epd(args.epd_bed_file, args.epd_motif_file)

    # Load 2nd EPD files if specified
    if((args.epd_bed_file_2 is not None) & (args.epd_motif_file_2 is not None) ):
      epd_motif = pd.concat([epd_motif_1, load_epd(args.epd_bed_file_2, args.epd_motif_file_2)])
    else:
      epd_motif = epd_motif_1

    # Load genes' information
    df_gene = get_gene_info(args.genes)

    # Pyranges for examining overlaps
    gr_promoter = df_to_pyranges(df_gene, start_col='tss', end_col='tss', start_slop=1000, end_slop=1000, strand_col='strand')
    gr_epd = df_to_pyranges(epd_motif, strand_col='strand')

    # Combine the promoter regions with EPD promoters
    df_combined = gr_promoter.join(gr_epd, strandedness='same', how='left').df[['chr', 'Start', 'End', 'strand', 'gene_name', 'gene_body_length', 'epd_id', 'TATA-box', 'Inr']]
    df_combined.replace({'epd_id': {'-1': ''}, 'TATA-box': {-1: 0}, 'Inr': {-1: 0}}, inplace=True)

    # Calculate burst sizes
    tata = df_combined['TATA-box']
    inr = df_combined['Inr']
    gl = np.log10(df_combined['gene_body_length'])
    # The equation and coefficients are from Larsson et al. (2019) Nature, 566, 251-254.
    df_combined['log10_burst_size'] = 0.865 - 0.0751 * gl + 0.6018 * tata - 0.1006 * gl * tata + 0.1516 * tata * inr

    # Output burst size information
    df_combined.to_csv(os.path.join(args.outdir, "epd_burst_sizes.txt"), sep="\t", index=False, header=True, float_format="%.6f")

    # Output gene-wise burst sizes
    df_out = df_combined[['gene_name', 'log10_burst_size']].groupby('gene_name', as_index=False)\
                        .aggregate({'log10_burst_size':['mean']})\
                        .set_axis(['gene_name', 'log10_burst_size'], axis='columns')
    df_out.to_csv(os.path.join(args.outdir, "gene_burst_sizes.txt"), sep="\t", index=False, header=False, float_format="%.6f")


if __name__ == '__main__':
    args = parseargs()
    main(args)


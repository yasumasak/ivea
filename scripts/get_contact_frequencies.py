import argparse
import pandas as pd
import numpy as np
import sys, traceback, os, os.path
import time
from tools import *
from hic import *

def parseargs():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description="Get enhancer-gene contact frequencies.", formatter_class=formatter)
    readable = argparse.FileType('r')

    # Basic parameters
    parser.add_argument('--enhancers', required=True, help="File of enhancer elements formatted as enhancer_elements.txt produced by annotate_regulatory_elements.py")
    parser.add_argument('--promoters', required=True, help="File of promoter elements formatted as gene_promoter_elements.txt produced by annotate_regulatory_elements.py")
    parser.add_argument('--outdir', required=True, help="Output directory")
    parser.add_argument('--window', type=int, default=3000000, help="Make contacts for all candidate elements within this distance of the gene's TSS")

    # HiC
    parser.add_argument('--hicdir', default=None, help="HiC directory: [hicdir]/chr*/chr*.*.gz. If not specified, computing contact frequencies from powerlaw only")
    parser.add_argument('--hic_resolution', type=int, help="HiC resolution")
    parser.add_argument('--tss_hic_contribution', type=float, default=100, help="Weighting of diagonal bin of hic matrix as a percentage of the maximum of its neighboring bins")
    parser.add_argument('--hic_pseudocount_distance', type=int, default=1e6, help="A pseudocount is added equal to the powerlaw fit at this distance")
    parser.add_argument('--hic_type', default="juicebox", choices=["juicebox","bedpe"], help="Format of hic files")

    # Power law
    parser.add_argument('--scale_hic_using_powerlaw', action="store_true", help="Scale Hi-C values using powerlaw relationship")
    parser.add_argument('--hic_gamma', type=float, default=.87, help="Powerlaw exponent of hic data. Must be positive")
    parser.add_argument('--hic_gamma_reference', type=float, default=.87, help="Powerlaw exponent to scale to. Must be positive")

    # Other
    parser.add_argument('--tss_slop', type=int, default=1000, help="Distance from tss to search for self-promoters")
    parser.add_argument('--chromosomes', default="all", help="Chromosomes to make contacts for. Defaults to intersection of all chromosomes in --promoters and --enhancers")
    parser.add_argument('--include_chrY', '-y', action="store_true", help="Make contacts on Y chromosome")

    args = parser.parse_args()
    return(args)


def main(args):

    print("Get enhancer-gene contact frequencies...")

    validate_args(args)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print_params(args)
    print("\n")

    print("  Reading promoter elements")
    promoters = pd.read_csv(args.promoters, sep = "\t")
    genes = promoters.loc[:,['chr','gene_name','tss']]
    genes.columns = ['chr','target_gene', 'target_tss']

    print("  Reading enhancer elements")
    enhancers_full = pd.read_csv(args.enhancers, sep = "\t")
    enhancers = enhancers_full.loc[:,['chr','start','end','name','class']]

    # Split chromosome list
    if args.chromosomes == "all":
        chromosomes = set(genes['chr']).intersection(set(enhancers['chr'])) 
        if not args.include_chrY:
            chromosomes.discard('chrY')
    else:
        chromosomes = args.chromosomes.split(",")

    # Make contacts for each chromosome
    for chromosome in chromosomes:
        print("  Making contacts for chromosome: {}".format(chromosome))
        t = time.time()

        this_enh = enhancers.loc[enhancers['chr'] == chromosome, :].copy()
        this_genes = genes.loc[genes['chr'] == chromosome, :].copy()

        contact = make_contacts(chromosome, this_enh, this_genes, args)

        # Output file
        file_contact = os.path.join(args.outdir, "enhancer-gene_contacts." + chromosome + ".txt")
        contact.to_csv(file_contact, sep="\t", index=False, header=True, float_format="%.6f")

        print("  Completed chromosome: {}. Elapsed time: {} \n".format(chromosome, time.time() - t))

    print("Done.\n")


def validate_args(args):
    if args.hicdir and args.hic_type == "juicebox":
        assert args.hic_resolution is not None, "HiC resolution must be provided if hic_type is juicebox"

    if not args.hicdir:
        print("WARNING: Hi-C directory is not provided. Computing contact frequencies from powerlaw only.")


def make_contacts(chromosome, enhancers, genes, args):
    cont = make_contact_table(chromosome, enhancers, genes, args)
    cont = annotate_contacts(cont, args.tss_slop)
    cont = add_powerlaw(cont, args)

    # If Hi-C directory is not provided, only powerlaw model will be computed
    if args.hicdir:
        hic_file, hic_norm_file, hic_is_vc = get_hic_file(chromosome, args.hicdir, hic_type = args.hic_type)
        cont = add_hic_to_enh_gene_table(enhancers, genes, cont, hic_file, hic_norm_file, hic_is_vc, chromosome, args)
        cont_cols = ['chr','start','end','name','target_gene','target_tss','distance','self_promoter','powerlaw_contact','powerlaw_contact_reference','hic_contact','hic_contact_scaled','hic_pseudocount','hic_contact_scaled_adjusted']
    else:
        cont_cols = ['chr','start','end','name','target_gene','target_tss','distance','self_promoter','powerlaw_contact','powerlaw_contact_reference']

    return cont[cont_cols]


def make_contact_table(chromosome, enh, genes, args):
    enh['enh_midpoint'] = (enh['start'] + enh['end'])/2
    enh['enh_idx'] = enh.index
    genes['gene_idx'] = genes.index
    enh_pr = df_to_pyranges(enh)
    genes_pr = df_to_pyranges(genes, start_col = 'target_tss', end_col = 'target_tss', start_slop=args.window, end_slop = args.window)

    cont = enh_pr.join(genes_pr).df.drop(['Start_b','End_b','chr_b','Chromosome','Start','End'], axis = 1)
    cont['distance'] = abs(cont['enh_midpoint'] - cont['target_tss'])
    cont = cont.loc[cont['distance'] < args.window,:] #for backwards compatability

    return cont


def annotate_contacts(cont, tss_slop=1000):
    cont['self_promoter'] = np.logical_and.reduce((cont['class'] == 'promoter' , cont.start - tss_slop < cont.target_tss, cont.end + tss_slop > cont.target_tss))

    return(cont)


def add_powerlaw(cont, args):
    cont['powerlaw_contact'] = get_powerlaw_at_distance(cont['distance'].values, args.hic_gamma)
    cont['powerlaw_contact_reference'] = get_powerlaw_at_distance(cont['distance'].values, args.hic_gamma_reference)

    return cont


def add_hic_to_enh_gene_table(enh, genes, cont, hic_file, hic_norm_file, hic_is_vc, chromosome, args):

    HiC = load_hic(hic_file = hic_file,
                    hic_norm_file = hic_norm_file,
                    hic_is_vc = hic_is_vc,
                    hic_type = args.hic_type,
                    hic_resolution = args.hic_resolution,
                    tss_hic_contribution = args.tss_hic_contribution,
                    window = args.window,
                    min_window = 0,
                    gamma = args.hic_gamma)

    # Add HiC data to 'cont' table
    # At this point we have a table where each row is an enhancer/gene pair.
    # We need to add the corresponding HiC matrix entry.
    # If the HiC is provided in juicebox format (ie constant resolution), then we can just merge using the indices
    # But more generally we do not want to assume constant resolution. In this case HiC should be provided in bedpe format

    if args.hic_type == "bedpe":
        # Use pyranges to compute overlaps between enhancers/genes and hic bedpe table
        # Consider each range of the hic matrix separately - and merge each range into both enhancers and genes. 
        # Then remerge on hic index

        HiC['hic_idx'] = HiC.index
        hic1 = df_to_pyranges(HiC, start_col='x1', end_col='x2', chr_col='chr1')
        hic2 = df_to_pyranges(HiC, start_col='y1', end_col='y2', chr_col='chr2')

        # Overlap in one direction
        enh_hic1 = df_to_pyranges(enh, start_col = 'enh_midpoint', end_col = 'enh_midpoint', end_slop = 1).join(hic1).df
        genes_hic2 = df_to_pyranges(genes, start_col = 'target_tss', end_col = 'target_tss', end_slop = 1).join(hic2).df
        ovl12 = enh_hic1[['enh_idx','hic_idx','hic_contact']].merge(genes_hic2[['gene_idx', 'hic_idx']], on = 'hic_idx')

        # Overlap in the other direction
        enh_hic2 = df_to_pyranges(enh, start_col = 'enh_midpoint', end_col = 'enh_midpoint', end_slop = 1).join(hic2).df
        genes_hic1 = df_to_pyranges(genes, start_col = 'target_tss', end_col = 'target_tss', end_slop = 1).join(hic1).df
        ovl21 = enh_hic2[['enh_idx','hic_idx','hic_contact']].merge(genes_hic1[['gene_idx', 'hic_idx']], on = ['hic_idx'])

        # Concatenate both directions and merge into preditions
        ovl = pd.concat([ovl12, ovl21]).drop_duplicates()
        cont = cont.merge(ovl, on = ['enh_idx', 'gene_idx'], how = 'left')
        cont.fillna(value={'hic_contact' : 0}, inplace=True)

    elif args.hic_type == "juicebox":
        # Merge directly using indices
        # Could also do this by indexing into the sparse matrix (instead of merge) but this seems to be slower
        # Index into sparse matrix
        # cont['hic_contact'] = [HiC[i,j] for (i,j) in cont[['enh_bin','tss_bin']].values.tolist()]

        cont['enh_bin'] = np.floor(cont['enh_midpoint'] / args.hic_resolution).astype(int)
        cont['tss_bin'] = np.floor(cont['target_tss'] / args.hic_resolution).astype(int)
        if not hic_is_vc:
            #in this case the matrix is upper triangular.
            #
            cont['bin1'] = np.amin(cont[['enh_bin', 'tss_bin']], axis = 1)
            cont['bin2'] = np.amax(cont[['enh_bin', 'tss_bin']], axis = 1)
            cont = cont.merge(HiC, how = 'left', on = ['bin1','bin2'])
            cont.fillna(value={'hic_contact' : 0}, inplace=True)
        else:
            # The matrix is not triangular, its full
            # For VC assume genes correspond to rows and columns to enhancers
            cont = cont.merge(HiC, how = 'left', left_on = ['tss_bin','enh_bin'], right_on=['bin1','bin2'])

        cont.fillna(value={'hic_contact' : 0}, inplace=True)

    cont.drop(['x1','x2','y1','y2','bin1','bin2','enh_idx','gene_idx','hic_idx','enh_midpoint','tss_bin','enh_bin'], inplace=True, axis = 1, errors='ignore')

    # Add powerlaw scaling
    cont = scale_hic_with_powerlaw(cont, args)

    # Add pseudocount
    cont = add_hic_pseudocount(cont, args)

    return(cont)


def scale_hic_with_powerlaw(cont, args):
    # Scale hic values to reference powerlaw

    if not args.scale_hic_using_powerlaw:
        cont['hic_contact_scaled'] = cont['hic_contact']
    else:
        cont['hic_contact_scaled'] = cont['hic_contact'] * (cont['powerlaw_contact_reference'] / cont['powerlaw_contact'])

    return(cont)


def add_hic_pseudocount(cont, args):
    # Add a pseudocount based on the powerlaw expected count at a given distance

    powerlaw_fit = get_powerlaw_at_distance(cont['distance'].values, args.hic_gamma)
    powerlaw_fit_at_ref = get_powerlaw_at_distance(args.hic_pseudocount_distance, args.hic_gamma)

    pseudocount = np.amin(pd.DataFrame({'a' : powerlaw_fit, 'b' : powerlaw_fit_at_ref}), axis = 1)
    cont['hic_pseudocount'] = pseudocount
    cont['hic_contact_scaled_adjusted'] = cont['hic_contact_scaled'] + pseudocount

    return(cont)



if __name__ == '__main__':
    args = parseargs()
    main(args)


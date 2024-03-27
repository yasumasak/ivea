import os
import numpy as np
import pandas as pd
from scipy import interpolate
import re
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError
import sys
import pyranges as pr
import linecache
import traceback
import time

pd.options.display.max_colwidth = 10000 #seems to be necessary for pandas to read long file names...

# setting this to raise makes sure that any dangerous assignments to pandas
# dataframe slices/subsets error rather than warn
pd.set_option('mode.chained_assignment', 'raise')

def run_command(command, **args):
    print("  Running command: " + command)
    return check_call(command, shell=True, **args)


def print_params(args):
    for arg in vars(args):
        print("  " + arg + " = " + str(getattr(args, arg)))


def df_to_pyranges(df, start_col='start', end_col='end', chr_col='chr', strand_col=None, start_slop=0, end_slop=0):
    df_pr = df.copy()
    df_pr['Chromosome'] = df[chr_col]
    df_pr['Start'] = df[start_col] - start_slop
    df_pr['Start'] = df_pr['Start'].where(df_pr['Start'] >=0, 0)
    df_pr['End'] = df[end_col] + end_slop
    if(strand_col is not None):
        df_pr['Strand'] = df[strand_col]

    return(pr.PyRanges(df_pr))


def sort_bed(bed_file, sorted_bed_file, genome_sizes, verbose=True):

    command = "bedtools sort -faidx {genome_sizes} -i {bed_file} > {sorted_bed_file}".format(**locals())

    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    if verbose: print("  Running: " + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')
    if not err == '':
        raise RuntimeError("Failed to sort {}.".format(bed_file))

    return stdoutdata


def run_count_reads(target, output, bed_file, genome_sizes, use_fast_count):
    if target.endswith(".bam"):
        count_bam(target, bed_file, output, genome_sizes=genome_sizes, use_fast_count=use_fast_count)
    else:
        raise ValueError("File {} name was not *.bam".format(target))


def count_bam(bamfile, bed_file, output, genome_sizes, use_fast_count=True, verbose=True):
    completed = True

    # Fast count:
    # bamtobed uses a lot of memory. Instead reorder bed file to match ordering of bam file. Assumed .bam file is sorted in the chromosome order defined by its header.
    # Then use bedtools coverage, then sort back to expected order
    # Requires an faidx file with chr in the same order as the bam file.
    if use_fast_count:
        temp_output = output + ".temp_sort_order"
        faidx_command = "awk 'FNR==NR {{x2[$1] = $0; next}} $1 in x2 {{print x2[$1]}}' {genome_sizes} <(samtools view -H {bamfile} | grep SQ | cut -f 2 | cut -c 4- )  > {temp_output}".format(**locals())
        command = "bedtools sort -faidx {temp_output} -i {bed_file} | bedtools coverage -g {temp_output} -counts -sorted -a stdin -b {bamfile} | awk '{{print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $NF}}'  | bedtools sort -faidx {genome_sizes} -i stdin > {output}; rm {temp_output}".format(**locals()) #

        #executable='/bin/bash' needed to parse < redirect in faidx_command
        p = Popen(faidx_command, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
        if verbose: print("  Running: " + faidx_command)
        (stdoutdata, stderrdata) = p.communicate()
        err = str(stderrdata, 'utf-8')

        p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        if verbose: print("  Running: " + command)
        (stdoutdata, stderrdata) = p.communicate()
        err = str(stderrdata, 'utf-8')

        try:
            data = pd.read_table(output, header=None).loc[:,3].values
        except Exception as e:
            print("Fast count method failed to count: " + str(bamfile) + "\n")
            print(err)
            print("Trying bamtobed method ...\n")
            completed = False

    # Alternate counting method. Slower and requires more memory.
    # convert BAM to BED, filter to standard chromosomes, sort, then use the very fast bedtools coverage -sorted algorithm
    # Note: This requires that bed_file is also sorted and in same chromosome order as genome_sizes (first do bedtools sort -i bed_file -faidx genome_sizes)
    #         BEDTools will error out if files are not properly sorted
    # Also requires that {genome_sizes} has a corresponding {genome_sizes}.bed file
    if not use_fast_count or ("terminated"  in err) or ("Error" in err) or ("ERROR" in err) or not completed:
        command = "bedtools bamtobed -i {bamfile} | cut -f 1-3 | bedtools intersect -wa -a stdin -b {genome_sizes}.bed | bedtools sort -i stdin -faidx {genome_sizes} | bedtools coverage -g {genome_sizes} -counts -sorted -a {bed_file} -b stdin | awk '{{print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $NF}}' > {output}".format(**locals())
        p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        if verbose: print("  Running: " + command)
        (stdoutdata, stderrdata) = p.communicate()

        try:
            data = pd.read_table(output, header=None).loc[:,3].values
        except Exception as e:
            print(e)
            print(stderrdata)
            completed = False

    # Check for successful finish -- BEDTools can run into memory problems
    err = str(stderrdata, 'utf-8')
    if ("terminated" not in err) and ("Error" not in err) and ("ERROR" not in err) and any(data):
        print("  Counting completed successfully. \n")
        completed = True
    else:
        print("Failed to count file: " + str(bamfile) + "\n")
        print(err)
        completed = False


def add_counts_for_bed(df_orig, bed_file, genome_sizes, feature_bam, directory, filebase, force=False, use_fast_count=True):
    orig_shape = df_orig.shape[0]
    feature_name = os.path.basename(feature_bam)
    counts_bedgraph = os.path.join(directory, "{}.{}.read_counts.bedgraph".format(filebase, feature_name))

    if force or (not os.path.exists(counts_bedgraph)) or (os.path.getsize(counts_bedgraph) == 0):
        print("  Regenerating", counts_bedgraph)
        print("  Counting coverage for {}".format(filebase + "." + feature_name))
        run_count_reads(feature_bam, counts_bedgraph, bed_file, genome_sizes, use_fast_count)
    else:
        print("  Loading coverage from pre-calculated file for {}".format(filebase + "." + feature_name))

    df_counts = read_bed6(counts_bedgraph)
    score_column = df_counts.columns[-1]

    df_counts = df_counts[['chr', 'start', 'end', score_column]]
    df_counts.rename(columns={score_column: "read_count"}, inplace=True)
    df_counts['chr'] = df_counts['chr'].astype('str')

    df_merged = df_orig.merge(df_counts.drop_duplicates())

    assert df_merged.shape[0] == orig_shape, "Dimension mismatch"

    return df_merged


def read_bed6(filename, chr=None, sort=False, skip_chr_sorting=True):
    skip = 1 if ("track" in open(filename, "r").readline()) else 0
    col_names = ["chr", "start", "end", "name", "score", "strand"]
    df = pd.read_table(filename, header=None, skiprows=skip, comment='#')
    n_cols = len(df.columns)
    if n_cols > 6:
        df = df.iloc[:, 0:6]
        df.columns = col_names
    else:
        df.columns = col_names[0:n_cols]
        
    df = df.dropna(axis=1, how='all')  # drop empty columns
    assert df.columns[0] == "chr"

    df['chr'] = pd.Categorical(df['chr'], ordered=True)
    if chr is not None:
        df = df[df.chr == chr]
    if not skip_chr_sorting:
        df.sort_values("chr", inplace=True)
    if sort:
        df.sort_values(["chr", "start", "end"], inplace=True)
    return df


def get_tss_for_bed(bed):
    assert_bed6(bed)
    tss = bed['start'].copy()
    tss.loc[bed.loc[:,'strand'] == "-"] = bed.loc[bed.loc[:,'strand'] == "-",'end'] - 1
    return tss


def assert_bed6(df):
    assert(type(df).__name__ == "DataFrame")
    assert('chr' in df.columns)
    assert('start' in df.columns)
    assert('end' in df.columns)
    assert('strand' in df.columns)


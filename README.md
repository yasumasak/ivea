
# IVEA : Integrative variational Bayesian inference of regulatory element activity for predicting enhancer–gene regulatory interactions

The IVEA predicts enhancer-gene regulatory interactions by estimating promoter and enhancer activities. This repository contains the scripts to prepare necessary input files as well as the scripts for the inference, and example commands for K562 cell line data set.


## Outline of workflow

The following data are used as input: 
* chromatin accessibility (DNase-seq or ATAC-seq)
* chromatin contact frequency (Hi-C) (Optional)
* gene expression (RNA-seq)
* genomic-sequence based burst sizes (Optional)

The scripts below in ```${IVEA_HOME}/scripts/``` proccess the above data into appropriate data format.
(```${IVEA_HOME}``` denotes a home directory of IVEA.)
* ```get_regulatory_elements.py```
* ```get_contact_frequencies.py```
* ```map_gene_expressions.py```
* ```get_burst_sizes.py```

  chromosome and gene annotation files (provided in ```${IVEA_HOME}/reference/``` for hg19) are used in the scripts.

Finally, ```run_IVEA.R``` in ```${IVEA_HOME}/scripts/``` performs variational Bayesian inference to predict enhancer-gene regulatory interactions.

The core functions for the variational inference are implemented in R scripts in ```${IVEA_HOME}/R``` directory as IVEA package. 
Users need to install the IVEA package (```R CMD INSTALL --no-multiarch --with-keep.source ${IVEA_HOME}```) before running the R script.

The example commands in the following sections are supposed to run in ```${IVEA_HOME}/example/``` and available in ```${IVEA_HOME}/example/commands_example.sh```


## Dependencies

The codebase relies on the following dependancies (tested version provided in parentheses):

```
Python (3.9.5)
R (3.5.0)
samtools (1.9)
bedtools (2.29.1)
MACS2 (2.2.7.1) - Partial dependancy
RSEM (RSEM-1.3.0_STAR-2.5.3a) - Partial dependancy
liftOver (2006-04-26) - Partial dependancy

Python packages:
pandas (1.2.5)
numpy (1.21.1)
scipy (1.7.0)
pyranges (0.0.113)
gtfparse (1.2.1)

R packages:
optparse (1.7.1)
data.table (1.14.0)
Matrix (1.2.14)
ghyp (1.6.1)
invgamma (1.1)
methods (3.5.0)
utils (3.5.0)
```


## Regulatory element

Regulatory elements (REs) are defined by peaks on a DNase-seq or ATAC-seq. The counts of reads mapped on the defined REs and the lengthes of REs are used in the variational inference.

Here we first use MACS2 to call peaks and then use ```get_regulatory_elements.py``` that defines REs, count mapped reads, and classify them into promoter/enhancer.

### Call peaks with MACS2

Example with K562 chr22:
```
macs2 callpeak \
-t ./input/wgEncodeUwDnaseK562AlnMerged.chr22.sorted.bam \
-n wgEncodeUwDnaseK562AlnMerged.chr22.macs2 \
-f BAM -g hs -p .1 \
--call-summits \
--outdir ./output
```

The resultant ```${dnase_name}_peaks.narrowPeak``` is used in the next.

### Define regulatory elements

The following processing steps are taken in ```get_regulatory_elements.py```:
 1. Count DNase-seq/Atac-seq reads in each peak and retain the top N peaks (```--n_peaks```) with the most read counts.
 2. Resize each of these N peaks to be a fixed number of base pairs (```--extension_from_summit```) centered on the peak summit.
 3. Remove any regions listed in the 'blacklist' (```--regions_blacklist```) and include any regions listed in the 'whitelist' (```--regions_whitelist```).
 4. Merge any overlapping regions. The merged peaks are defiend as REs.
 5. Classify the REs. A promoter element for a gene is defined as sum of REs in a region that spans a fixed number of base pairs (```--tss_slop_for_class_assignment```) from the gene TSSs. All REs are considered as enhancer elements, including REs found in the promoter regions, since gene promoters can potentially act as enhancers (Dao et al., 2017; Andersson and Sandelin, 2020).

Example with K562 chr22:
```
python ${IVEA_HOME}/scripts/get_regulatory_elements.py \
--outdir ./output \
--narrowPeak ./output/wgEncodeUwDnaseK562AlnMerged.chr22.macs2_peaks.narrowPeak \
--bam ./input/wgEncodeUwDnaseK562AlnMerged.chr22.sorted.bam \
--chrom_sizes ${IVEA_HOME}/reference/hg19.chrom.sizes \
--regions_blacklist ${IVEA_HOME}/reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
--extension_from_summit 250 \
--genes ${IVEA_HOME}/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.excl_chrY.Fulco_2019.bed \
--tss_slop_for_class_assignment 1000 \
--n_peaks 3000 # 150000 (default) is recommended for genome-wide analysis. 3000 is just used for the small example on chr22.
```
Main outputs:
* **enhancer_elements.txt**: Enhancer elements with Dnase-seq (or ATAC-seq) read counts.
* **gene_promoter_elements.txt**: Gene promoter elements with Dnase-seq (or ATAC-seq) read counts.


## Contact frequency

Contact frequencies between the gene TSSs and the enhancer elements are used in the variational inference. ```get_contact_frequencies.py``` processes Hi-C data with gene TSSs and enhancer elements information to provide contact frequencies of enhancer-gene pairs.

The following Hi-C data processing steps reported in Fulco et al (2019) are used:
1. Each diagonal entry of the Hi-C matrix is replaced by the maximum of its four neighbouring etries.
2. All entries of the Hi-C matrix with a value of NaN or corresponding to KR normalization factors < 0.25 are replaced with the expected contact under the power-law distribution with the law's exponent (```--hic_gamma```).
3. A small adjustment (pseudocount) is added to the entries of the Hi-C matrix. For the entries with distance larger than the pseudocount distance (```--hic_pseudocount_distance```), the expected contact frequency under the power-law distribution is added. For those within the pseudocount distance, a constant adjustment equal to the expected contact frequency at the pseudocount distance is added.

### Format of Hi-C data
* Juicer format: Three column 'sparse matrix' format representation of a Hi-C matrix.
* BEDPE format: More general format which can support variable and arbitrary bin sizes by specifying ```--hic_type bedpe```. Note that if contact data is provided in BEDPE format, the 1st and 2nd processing of the Hi-C data described above are **not** applied. The BEDPE file should be a tab-delimited file containing 8 columns (chr1,start1,end1,chr2,start2,end2,name,score) where score denotes the contact frequency. 

### Without experimental Hi-C contact data
If experimentally derived contact data is not available, two alternative approaches can be taken.
* Powerlaw-estimate: The powerlaw-estimated contact frequency is applied by not specifying ```--hicdir```. It has been shown that Hi-C contact frequencies generally follow a powerlaw relationship (with respect to genomic distance) and that many TADs, loops and other structural features of the 3D genome are **not** cell-type specific (Sanborn et al 2015, Rao et al 2014). 
* Average Hi-C: The average Hi-C matrix (averaged across 10 cell lines, at 5kb resolution: GM12878, NHEK, HMEC, RPE1, THP1, IMR90, HUVEC, HCT116, K562, KBM7) can be downloaded from: <ftp://ftp.broadinstitute.org/outgoing/lincRNA/average_hic/average_hic.v2.191020.tar.gz> (20 GB). The average Hi-C profile showed approximately equally good performance as using a cell-type specific Hi-C profile (Fulco et al 2019). 

Example with K562 chr22 without Hi-C contact data:
```
python ${IVEA_HOME}/scripts/get_contact_frequencies.py \
--enhancers ./output/enhancer_elements.txt \
--promoters ./output/gene_promoter_elements.txt \
--window 5000000 \
--outdir ./output \
--chromosomes chr22 \
#--hicdir ./input/HiC/raw \ # Set when using Hi-C contact data.
#--hic_resolution 5000 \ # Set when using Hi-C contact data.
```
Main outputs:
* **enhancer-gene_contacts.chr?.txt**: Contact frequencies between the gene TSSs and the enhancer elements.


## Gene expression

Gene-level RNA-seq read counts and effective lengths are used in the variational inference. The genes analyzed in the variational inference are filtered based on transcripts per million (TPM) by default.

### RSEM 

Here we use ```'genes.result'``` generated by RSEM that contains gene-level read counts, TPM and effective lengths.

Example with K562:
```
rsem-calculate-expression
--star-gzipped-read-file --no-bam-output --star-output-genome-bam --estimate-rspd
--star --star-path STAR_PATH
--paired-end ./input/ENCFF001REG.fastq.gz ./input/ENCFF001REF.fastq.gz
${reference_gencode_v26lift37} ./output/ENCFF001REG-ENCFF001REF_rsem
```
```${reference_gencode_v26lift37}``` is a Gencode-based reference generated by RSEM ```rsem-prepare-reference``` command using ```${IVEA_HOME}/reference/gencode.v26lift37.annotation.gtf```.

### Refseq-based reference

In the K562 example, we use RefSeq-based gene annotation as a reference which is different from one used in RSEM (Gencode-based). In such case, ```map_gene_expressions.py``` can be used to map Gencode-based RSEM ```'genes.result'``` to the RefSeq-based reference. 

Example with K562:
```
python ${IVEA_HOME}/scripts/map_gene_expressions.py \
--bed_ref ${IVEA_HOME}/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.excl_chrY.Fulco_2019.bed \
--gtf_gencode ${IVEA_HOME}/reference/gencode.v26lift37.annotation.gtf.gz \
--rsem ./output/ENCFF001REG-ENCFF001REF_rsem.genes.results \
--outdir ./output
```
Main outputs:
* **gene_expressions.txt**: Gene-level read counts, TPM and effective lengths.

### Gencode-based reference

In case of using Gencode-based reference (same as in RSEM) in IVEA, ```get_gencode_bed.py``` can be used to generate gene position BED file and gene id/name list file from the Gencode gtf file, and ```get_gencode_expression.py``` can be used to make ```'gene_expressions.txt'``` from RSEM ```'genes.result'```. 


## Transcriptional burst size

Transcriptional burst size estimate is used to scale the gene promoter activity in the variational inference. Larsson et al (2019) found that burst size can be estimated from core promoter sequence elements and gene body length. ```get_burst_sizes.py``` utilizes the regression formula reported in Larsson et al (2019) and provides gene-wise burst size estimates. 

Firstly, the Eukaryotic Promoter Database (EPD) data are needed to be downloaded from ftp://ccg.epfl.ch/.
For the human reference genome hg19, the following files were downloaded to ```${EPD_dir}```.
- epdnew/H_sapiens/005/Hs_EPDnew_005_hg19.bed
- epdnew/H_sapiens/005/db/promoter_motifs.txt
- epdnew/H_sapiens_nc/001/HsNC_EPDnew_001_hg38.bed
- epdnew/H_sapiens_nc/001/db/promoter_motifs.txt
The HsNC_EPDnew_001_hg38.bed was lifted to hg19 by liftOver.

Example for hg19:
```
python ${IVEA_HOME}/scripts/get_burst_sizes.py \
--epd_bed_file ${EPD_dir}/epdnew/H_sapiens/005/Hs_EPDnew_005_hg19.bed \
--epd_motif_file ${EPD_dir}/epdnew/H_sapiens/005/db/promoter_motifs.txt \
--epd_bed_file_2 ${EPD_dir}/epdnew/H_sapiens_nc/001/HsNC_EPDnew_001_hg19.lifted.bed \
--epd_motif_file_2 ${EPD_dir}/epdnew/H_sapiens_nc/001/db/promoter_motifs.txt \
--genes ${IVEA_HOME}/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.excl_chrY.Fulco_2019.bed \
--outdir ./output
```
Main outputs:
* **gene_burst_sizes.txt**: Gene-wise burst size estimates.

A set of burst sizes obtained by the example command above for hg19 is available as ```${IVEA_HOME}/reference/gene_burst_sizes.hg19.txt```. 


## Variational Bayesian inference for predicting enhancer-gene regulatory interactions

The variational Bayesian inference for predicting enhancer-gene regulatory interactions is made by  ```run_IVEA.R```. It provides estimates of promoter and enhancer activities, and scores of enhancer-gene regulatory interactions.

The core functions for the variational inference are implemented in R scripts in ```${IVEA_HOME}/R``` directory as IVEA package. 
Users need to install the IVEA package (```R CMD INSTALL --no-multiarch --with-keep.source ${IVEA_HOME}```) before running the R script.

The inference can be run in parallel in a chromosomal basis. 

Example with K562 chr22:
```
Rscript ${IVEA_HOME}/scripts/run_IVEA.R \
--chr_region chr22 \
--alpha_pro 80 \
--alpha_enh 10 \
--contacts ./output/enhancer-gene_contacts.chr22.txt \
--rnas ./output/gene_expressions.txt \
--rna_cutoff_value 8 \
--enhancers ./output/enhancer_elements.txt \
--promoters ./output/gene_promoter_elements.txt \
--burst_sizes ${IVEA_HOME}/reference/gene_burst_sizes.hg19.txt \
--outdir ./output
```
Main outputs:
* **predictions_score.chr?.bedpe**: Prediction result in BEDPE format (enhancer and gene TSS position, and their interaction score).
* **predictions_info.chr?.txt**: Prediction result with detailed information: gene (name, chromosome, TSS), promoter (name, read count, length, activity), enhancer (name, read count, length, activity), and regulatory interaction (distance, contact frequency, strength, contribution, score).
* **estimates.enhancer_activity.chr?.bed**: Estimated enhancer activities in BED format (enhancer position, name and estimated enhancer activity).
* **estimates.promoter_activity.chr?.bed**: Estimated promoter activities in BED format (gene TSS, name and estimated promoter activity).


## Contact
Please submit a github issue with any questions or if you experience any issues/bugs. 

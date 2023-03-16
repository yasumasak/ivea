
IVEA_HOME=".."


## Regulatory element
### Define regulatory elements
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


## Contact frequency
python ${IVEA_HOME}/scripts/get_contact_frequencies.py \
--enhancers ./output/enhancer_elements.txt \
--promoters ./output/gene_promoter_elements.txt \
--window 3000000 \
--outdir ./output \
--chromosomes chr22 \
#--hicdir ./input/HiC/raw \ # Set when using Hi-C contact data.
#--hic_resolution 5000 \ # Set when using Hi-C contact data.


## Gene expression
### Refseq-based reference
python ${IVEA_HOME}/scripts/map_gene_expressions.py \
--bed_ref ${IVEA_HOME}/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.excl_chrY.Fulco_2019.bed \
--gtf_gencode ${IVEA_HOME}/reference/gencode.v26lift37.annotation.gtf.gz \
--rsem ./output/ENCFF001REG-ENCFF001REF_rsem.genes.results \
--outdir ./output


## Variational Bayesian inference for predicting enhancer-gene regulatory interactions
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


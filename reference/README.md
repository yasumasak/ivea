# IVEA reference data

IVEA itself is **genome-agnostic** — the R package (`R/*.R`), the driver
(`scripts/run_IVEA.R`), and the Python helpers carry no hardcoded build string,
chromosome sizes, or coordinate offsets. The genome build is determined entirely by
the reference files you pass on the command line. The only real assumption is UCSC
`chr`-prefixed chromosome names, which both hg19 and hg38 (UCSC) satisfy.

This directory ships two consistent reference sets:

- **hg19** — the build used for the analyses in Kimura et al. 2024. This is the
  reproduce-the-paper default and backs the bundled chr22 K562 example.
- **hg38** — an additive, batteries-included set for new studies aligned to GRCh38.

Files are named with a build suffix and live flat in this directory. Mixing files
from different builds in one run will silently produce wrong results — keep every
input on the same build.

## Build map

| File | Build | Consumer (CLI arg) | Source / provenance | How obtained |
|---|---|---|---|---|
| `hg19.chrom.sizes` (+`.bed`) | hg19 | `get_regulatory_elements.py` (`--chrom_sizes`) | UCSC | committed |
| `hg38.chrom.sizes` (+`.bed`) | hg38 | `get_regulatory_elements.py` (`--chrom_sizes`) | UCSC (25 primary contigs chr1–22,X,Y,M, size-descending) | committed |
| `wgEncodeHg19ConsensusSignalArtifactRegions.bed` | hg19 | `get_regulatory_elements.py` (`--regions_blacklist`) | ENCODE hg19 consensus signal artifact regions | committed |
| `hg38-blacklist.v2.bed` | hg38 | `get_regulatory_elements.py` (`--regions_blacklist`) | ENCODE Blacklist **v2** (Amemiya et al. 2019; `ENCFF356LFX`, Boyle-Lab/Blacklist), first 3 columns, sorted | committed |
| `RefSeqCurated.170308.bed.CollapsedGeneBounds.excl_chrY.Fulco_2019.bed` | hg19 | `get_regulatory_elements.py` (`--genes`), `map_gene_expressions.py` (`--bed_ref`), `get_burst_sizes.py` (`--genes`) | RefSeq-curated, collapsed to one interval per gene (Fulco et al. 2019 / ABC model) | committed |
| `CollapsedGeneBounds.hg38.bed` | hg38 | same as the hg19 collapsed gene bounds above | the hg19 file lifted to hg38 (UCSC `hg19ToHg38.over.chain`, then filtered to the 25 primary contigs; 24,396 genes, unique names, BED-6) | committed |
| `gencode.v26lift37.annotation.gtf.gz` | hg19 | RSEM reference; `map_gene_expressions.py` (`--gtf_gencode`); `get_gencode_bed.py` | GENCODE release 26, GRCh37-mapped | committed |
| `gencode.v26.annotation.gtf.gz` | hg38 | same as the hg19 GTF above | GENCODE release 26, native GRCh38 | **fetched** (`scripts/fetch_reference.sh hg38`; git-ignored) |
| `gene_burst_sizes.hg19.txt`, `epd_burst_sizes.hg19.txt` | hg19 | `run_IVEA.R` (`--burst_sizes`) | `get_burst_sizes.py`, EPD hg19 promoters | committed |
| `gene_burst_sizes.hg38.txt`, `epd_burst_sizes.hg38.txt` | hg38 | `run_IVEA.R` (`--burst_sizes`) | `get_burst_sizes.py`, EPD hg38 promoters (see below) | committed |

Burst sizes are **optional** to `run_IVEA.R` (it defaults each gene's burst size to 1
when `--burst_sizes` is omitted).

## Why the blacklist naming differs between builds

hg19 uses the older ENCODE "consensus signal artifact regions" product, which was
never produced for hg38. Its maintained successor covering hg38 is the ENCODE
Blacklist **v2** (Amemiya et al. 2019). The `.v2` is that product's canonical upstream
filename and is kept for provenance. hg19 is intentionally **not** switched to v2 —
doing so would change hg19 outputs and break paper reproduction.

## Fetched (non-committed) files

Only the GENCODE GTF is fetched on demand, because of its size:

```
bash scripts/fetch_reference.sh hg38   # -> reference/gencode.v26.annotation.gtf.gz
bash scripts/fetch_reference.sh hg19   # verifies the vendored hg19 GTF
```

### GENCODE release and why v26 is the default

Both builds default to GENCODE **release 26** — the release used in the paper. The
vendored hg19 GTF is `gencode.v26lift37` (release 26 mapped to GRCh37); the hg38
default is the native GRCh38 `gencode.v26`. `CollapsedGeneBounds.*.bed` and the
committed burst sizes are all built on this same release, so the provided reference
set is internally consistent by default.

`fetch_reference.sh` accepts an optional release argument to pull a newer GENCODE
release, e.g. `bash scripts/fetch_reference.sh hg38 48`. **If you use a release other
than 26, regenerate the gene-derived files from that same release so they match the
GTF's gene set**, otherwise a run mixes annotations:

- `CollapsedGeneBounds.hg38.bed` — rebuild from the new release (e.g. collapse the
  release's gene intervals per gene name), or lift a matching hg19 collapsed set.
- `gene_burst_sizes.hg38.txt` / `epd_burst_sizes.hg38.txt` — rerun `get_burst_sizes.py`
  with `--genes` pointing at the rebuilt gene bounds (EPD promoter files are
  release-independent and can be reused).

IVEA itself imposes no release requirement; consistency across the GTF, gene bounds,
RSEM reference, and burst sizes is what matters.

## Regenerating the hg38 burst sizes

`gene_burst_sizes.hg38.txt` / `epd_burst_sizes.hg38.txt` were produced from EPD hg38
promoters. Download the EPD files (served over HTTPS at `https://epd.expasy.org/ftp/`):

- coding: `epdnew/H_sapiens/006/Hs_EPDnew_006_hg38.bed` + `epdnew/H_sapiens/006/db/promoter_motifs.txt`
- non-coding: `epdnew/H_sapiens_nc/001/HsNC_EPDnew_001_hg38.bed` + `epdnew/H_sapiens_nc/001/db/promoter_motifs.txt`

The EPD `.bed` files are whitespace-delimited; convert them to tab-delimited before
use (e.g. `awk 'BEGIN{OFS="\t"}{$1=$1;print}'`). Then:

```
python scripts/get_burst_sizes.py \
  --epd_bed_file   Hs_EPDnew_006_hg38.tab.bed  --epd_motif_file   promoter_motifs.txt \
  --epd_bed_file_2 HsNC_EPDnew_001_hg38.tab.bed --epd_motif_file_2 promoter_motifs.nc.txt \
  --genes reference/CollapsedGeneBounds.hg38.bed --outdir <outdir>
```

Unlike the hg19 workflow (which lifts the non-coding EPD file *down* to hg19), the
non-coding EPD file is already native hg38, so no liftOver is needed.

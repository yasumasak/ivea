# IVEA 1.2.0

Maintenance and infrastructure release. The variational Bayes inference and its
numerical results are unchanged from 1.1.0 — the single change to the inference
code (below) is numerically equivalent. This release modernizes the toolchain and
adds documentation, continuous integration, and hg38 support.

## Toolchain and dependencies

* The pipeline now runs on **Python 3.12** and **R 4.5**, and the whole stack
  installs in one step from the new `environment.yml`
  (`mamba env create -f environment.yml`).
* **Removed the `ghyp` dependency.** The generalized inverse Gaussian expectations
  in `R/openness.R` are now computed in closed form with base R `besselK`; the
  numerical-integration fallback for the overflow/divergence cases is retained.
  Results are numerically equivalent to 1.1.0.
* Python preprocessing modernized: `pyranges` → **bioframe** (pandas-native,
  supports pandas 2 / NumPy 2); `gtfparse` dropped in favour of a small built-in
  GTF reader; peak calling **MACS2 → MACS3**; NumPy 2 fix (`np.Inf` → `np.inf`).

## Genome support

* Added an **hg38** reference set (chrom sizes, blacklist, collapsed gene bounds,
  burst sizes) alongside the existing hg19 data, plus `scripts/fetch_reference.sh`
  to download the one large reference that is not vendored (the GENCODE GTF). hg19
  behaviour is unchanged.

## Documentation and citation

* Added **`CITATION.cff`** and **`inst/CITATION`** so IVEA is directly citable:
  GitHub's "Cite this repository" button and `citation("IVEA")` both return
  Kimura et al., *Bioinformatics Advances* 2024;4(1):vbae118
  (<https://doi.org/10.1093/bioadv/vbae118>).
* Rewrote the **README** (quick start, workflow schematic, "When to use IVEA").
* Scaffolded a **pkgdown** documentation site (`_pkgdown.yml`).
* Added a **K562 chr22 tutorial** (`vignettes/articles/ivea-k562.Rmd`): it runs the
  real inference and walks TPST2's predicted regulatory landscape with plots and
  the bundled CRISPRi-validated reference pairs.

## Continuous integration

* **`r-check.yml`** — `R CMD check` plus the `testthat` suite on R 4.5.
* **`pipeline.yml`** — an end-to-end check that re-runs the full chr22 example
  (Python preprocessing + R inference) into a clean directory and verifies the
  outputs reproduce the committed fixtures (`scripts/compare_outputs.R`).

## Example data and packaging

* Regenerated the chr22 `example/output/` fixtures so they are self-consistent
  with the current minus-strand TSS convention, and added the estimated
  promoter/enhancer activity output tables.
* `DESCRIPTION`: declared `URL` and `BugReports`, moved the tutorial tooling
  (`knitr`, `rmarkdown`, `ggplot2`) to `Suggests`, set `Depends: R (>= 4.0)`, and
  replaced `RoxygenNote` with `Config/roxygen2/version`. `.Rbuildignore` trims the
  build tarball from ~114 MB to ~19 KB.


# IVEA 1.1.0

* Added the **IVEA_nolE** model variant — the enhancer-activity prior without
  relative enhancer length — selectable via the `no_le_enh` switch.
* Output now reports **95% credible intervals and standard deviations** for the
  estimated promoter and enhancer activities, and writes the per-element activity
  estimate tables.
* Added the `testthat` test suite and the chr22 example outputs.


# IVEA 1.0.0

* First public release, accompanying Kimura et al. (2024), *Bioinformatics
  Advances*.

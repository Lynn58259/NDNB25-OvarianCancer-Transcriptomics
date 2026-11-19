# NDNB25 Treatment and Transcriptomic Profiling in Ovarian Cancer Models

This repository provides the processed inputs, analysis code, and numerical results for RNA-seq analyses associated with an IJMS manuscript investigating HSP90β inhibition by NDNB25 and SREBP1-associated lipid metabolism in ovarian cancer models. The associated manuscript, **Selective HSP90β Inhibition Disrupts SREBP1-Associated Metabolic Signaling and Suppresses Ovarian Cancer Progression: A Pilot Study**, is in proofreading. This repository covers the cell-line RNA-seq analyses for Figures 5 and 6; it does not reproduce the separate TCGA, animal, or protein experiments. Generated figure files are local outputs and are not distributed in this repository.

## Experimental design

OVCAR5 and OVCAR8 human ovarian cancer cells and ID8-Trp53−/− mouse ovarian cancer cells were treated with vehicle or 400 nM NDNB25 for 24 hours. Each cell line has three independent biological replicates per treatment (18 samples in total). The analysis pairs vehicle and treatment samples within three experimental batches per cell line. The original submission sample list confirms the sample-number mapping, and the GEO submission workbook agrees with all 18 batch labels. Neither source documents passage or collection dates, so independent laboratory-record confirmation of the biological pairing remains required. These labels must not be interpreted as documented library-preparation or sequencing batches. Plasmidsaurus performed 3′-end, UMI-based RNA sequencing. Sample assignments and GEO library names are recorded in [`data/metadata/samples.csv`](data/metadata/samples.csv).

## Data availability and inputs

The study is associated with NCBI GEO Series [GSE349387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349387). Public access to the accession could not be independently confirmed when this README was prepared; verify its release status before making this repository public. GEO holds the submitted raw sequencing data. This repository begins with the processed provider expression matrix, [`data/counts/provider_expression_matrix.tsv`](data/counts/provider_expression_matrix.tsv), which retains the provider's fractional counts, CPM values, gene annotations, and original sample IDs. The workflow selects species-specific count columns, rounds counts to integers, and filters genes before DESeq2. It does not align FASTQ reads or perform UMI counting. See [`DATA_PROVENANCE.md`](DATA_PROVENANCE.md) for the source and gene-set details.

## Analysis

The workflow runs DESeq2 separately for OVCAR5, OVCAR8, and ID8. Its design, `~ batch + treatment`, accounts for the paired batches while estimating the NDNB25-versus-vehicle effect within each cell line. Variance-stabilized counts are used for PCA. Ranked differential-expression results feed Hallmark and Gene Ontology Biological Process (GO BP) GSEA. A prespecified SREBP1-associated lipogenic module and selected lipid-metabolism programs and genes support the metabolic signaling analysis. The code generates the analyses corresponding to Figures 5 and 6 of the associated manuscript without changing the validated statistical methods.

## Repository structure

```text
.
├── data/                 # provider matrix, metadata, frozen gene sets
├── scripts/              # preprocessing, DESeq2, PCA, GSEA, plotting
├── results/
│   ├── differential_expression/  # complete per-cell-line results and VST
│   ├── normalized_counts/        # per-cell-line normalized count matrices
│   ├── pathway_enrichment/       # Hallmark, GO BP, lipogenic module
│   └── plot_data/                # PCA and figure source tables
├── environment/          # dependency installer and recorded versions
├── run_analysis.R         # complete workflow entry point
├── DATA_PROVENANCE.md
├── CITATION.cff
└── LICENSE
```

## Reproducing the analysis

Run all commands **from the repository root**, in a fresh clone or disposable copy. The workflow overwrites numerical results under `results/`; preserve a separate reference checkout before running if you want to compare outputs. Committed inputs and reference environment records are not modified by the entry point.

The validated environment uses **R 4.5.1, Bioconductor 3.22, Python 3.12.14, NumPy 2.3.5, and pandas 3.0.1**. Direct R dependencies are data.table, DESeq2, fgsea, ggplot2, patchwork, ggrepel, ggnewscale, ragg, svglite, and scales. BiocManager installs Bioconductor dependencies. Python dependencies are pinned in [`environment/requirements.txt`](environment/requirements.txt). Recorded R versions and session information remain in `environment/`.

1. Install R and Python. Put `Rscript` and `python` on `PATH`. On Windows, disable or bypass a nonfunctional Windows Store Python alias. If using another executable, set `PYTHON` to its full path before running R (PowerShell: `$env:PYTHON = "C:/path/to/python.exe"`; POSIX shell: `export PYTHON=/path/to/python`).
2. Install Python dependencies into an isolated environment, then select that environment's Python executable:

   ```sh
   python -m venv .venv
   # Activate .venv using your shell's normal activation command, then:
   python -m pip install -r environment/requirements.txt
   ```

3. **Standard method (recommended).** renv is not required. The existing installer is a convenience path; it installs missing packages and pins the Bioconductor release, but **does not enforce every recorded package version**:

   ```sh
   Rscript --vanilla environment/install_packages.R
   Rscript --vanilla run_analysis.R
   ```

   **Optional renv method (clean restoration not validated).** Skip this method to use the standard workflow above. An exact-version specification is provided in [`renv.lock`](renv.lock), containing only the workflow dependencies and their required dependency closure. To attempt an isolated restore:

   ```sh
   Rscript --vanilla -e "install.packages('renv', repos='https://cloud.r-project.org')"
   Rscript --vanilla environment/restore_locked.R
   Rscript --vanilla -e ".libPaths(c('renv/library', .libPaths())); source('run_analysis.R')"
   ```

   **Clean-install limitation:** the November 18, 2026 audit attempted a cache-disabled restore into an empty library. It failed while installing compiled packages, and a compiler probe confirmed missing Windows build tools (`'make' not found`). Restore also emitted binary-repository index warnings; missing build tools have not been established as the sole obstacle. Installing the appropriate [Rtools for R 4.5](https://cran.r-project.org/bin/windows/Rtools/) and retesting is required before claiming successful clean installation. On Linux/macOS, source packages can require compilers and graphics/font system libraries. A lockfile records package versions; it does not install R, Python, compilers, fonts, or system libraries. This optional method is provided without a clean-install guarantee and is not a prerequisite for the standard workflow or public release. Neither the standard installer nor the analysis entry point loads renv. There is no automatic activation, startup hook, or required renv bootstrap file; the restore helper is invoked explicitly. Local libraries and caches are ignored by Git, while `renv.lock` and the helper are retained. Existing package records are retained.

4. The entry point executes scripts 01–06: preprocessing, DESeq2/VST, PCA, Hallmark/GO BP GSEA, the lipogenic module, and plotting. No network access or external experimental data files are used during analysis. Installations require network access. Numerical outputs are written under `results/`; PDF, SVG, PNG, and figure source CSV files are generated under the Git-ignored `figures/` directory. Runtime package/session records are written under `results/quality_control/`.
5. Compare a regenerated checkout with a separate untouched reference checkout:

   ```sh
   python scripts/compare_results.py --reference ../reference_checkout --generated . --output comparison.csv
   ```

   This compares every committed numerical CSV in `results/` (excluding local quality-control outputs), checks schema, row order, missing values, categorical cells, and significance thresholds, and uses absolute tolerance 1e-10 plus relative tolerance 1e-8 for numeric cells. Any mismatch returns a nonzero exit status. PDF/SVG metadata and graphics rendering are not expected to be byte-identical across systems.

Validation regenerated all 36 reference numerical CSVs with exact numerical agreement using the documented R/Python versions. Clean renv restoration remains unvalidated, as explained above. Internal audit reports are retained locally and are not distributed with this repository.

## Method details and limitations

- Counts are selected by species-specific Ensembl prefix (`ENSG` for human, `ENSMUSG` for mouse), rounded using R's ties-to-even `round()`, and retained when at least 10 rounded counts occur in at least three of the six samples. Rounding does not recover uniquely assigned molecules from fractional provider counts.
- Each cell line uses `~ batch + treatment`, with CON as reference and NB25 versus CON as contrast. DESeq2 uses its default size-factor normalization and Wald testing; BH-adjusted p < 0.05 defines differential expression. All 18 samples are retained. VST uses `blind=FALSE`; PCA centers the 1,000 most variable genes within each cell line without scaling to unit variance.
- GSEA ranks genes with nonmissing DESeq2 adjusted p values by Wald statistic, uppercases symbols, and retains the first occurrence of duplicate symbols. Hallmark/GO BP use sizes 10–500; the fixed module uses sizes 5–30. All fgsea calls use `eps=0`, `nproc=1`, and seeds 20260923–20260925 in OVCAR5/OVCAR8/ID8 order. PCA labels use a fixed ggrepel seed.
- Hallmark and GO BP BH corrections are separate within each cell line and collection. Module q values are BH-adjusted across the three prespecified cell-line tests. Selected Hallmark display pathways were chosen after exploration. These testing families must not be conflated.
- The SREBP1-associated module is an association-based transcriptional readout, not a direct binding or target-engagement assay. Human SCD maps to mouse Scd1 and Scd2; human Scd2 heatmap positions are structurally not applicable. The module contains seven human and eight projected mouse members.
- Frozen MSigDB 2026.1.Hs sets avoid the previously documented msigdbr mouse ortholog cache-order issue. The analysis does not install or call msigdbr or regenerate ortholog mappings. Small sample size, fractional-count rounding, and unverified independent batch provenance limit interpretation. The repository reproduces processed-matrix analyses, not FASTQ alignment or UMI counting.

## Finding the results

Complete DESeq2 tables, including nonsignificant genes, are in `results/differential_expression/<cell_line>/all_genes.csv`. Normalized counts are in `results/normalized_counts/<cell_line>/normalized_counts.csv`. Complete Hallmark and GO BP enrichment tables, including nonsignificant pathways, and SREBP1-associated module results are in `results/pathway_enrichment/<cell_line>/`. PCA coordinates, explained variance, selected-gene values, and other figure source tables are in `results/plot_data/`. Replace `<cell_line>` with `OVCAR5`, `OVCAR8`, or `ID8`.

## Citation and license

The associated IJMS manuscript is in proofreading. The author list was confirmed by the code owner; no publication DOI, publication year, volume, or pages are asserted. For the software repository, use [`CITATION.cff`](CITATION.cff), and cite GEO Series GSE349387 for the sequencing dataset when the record becomes public.

Original analysis code is released under the [MIT License](LICENSE). Third-party gene sets and other source materials retain their own licensing and attribution conditions, described in [`data/gene_sets/README.md`](data/gene_sets/README.md).

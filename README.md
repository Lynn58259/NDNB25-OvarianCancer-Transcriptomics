# NDNB25 Treatment and Transcriptomic Profiling in Ovarian Cancer Models

This repository provides the processed inputs, analysis code, and numerical results for RNA-seq analyses associated with an IJMS manuscript investigating HSP90β inhibition by NDNB25 and SREBP1-associated lipid metabolism in ovarian cancer models. The workflow reproduces the manuscript's transcriptomic profiling and metabolic signaling analyses. Generated figure files are local outputs and are not distributed in this repository.

## Experimental design

OVCAR5 and OVCAR8 human ovarian cancer cells and ID8-Trp53−/− mouse ovarian cancer cells were treated with vehicle or 400 nM NDNB25 for 24 hours. Each cell line has three independent biological replicates per treatment (18 samples in total). The analysis pairs vehicle and treatment samples within three batches per cell line; these batch assignments should be confirmed against the original experimental records. Plasmidsaurus performed 3′-end, UMI-based RNA sequencing. Sample assignments and GEO library names are recorded in [`data/metadata/samples.csv`](data/metadata/samples.csv).

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

Run the following commands **from this repository's root directory**. The validated run used R 4.5.1, Python 3.12.14, NumPy 2.3.5, and pandas 3.0.1. Python dependencies are pinned in [`environment/requirements.txt`](environment/requirements.txt); the R installer adds the required CRAN and Bioconductor packages. The recorded R package versions are in [`environment/package_versions.txt`](environment/package_versions.txt). R package installation needs network access; running the analysis uses the frozen local gene sets.

```sh
python -m pip install -r environment/requirements.txt
Rscript --vanilla environment/install_packages.R
Rscript --vanilla run_analysis.R
```

The entry point executes scripts 01–06 in order. It uses `python` on `PATH` for the PCA helper; set `PYTHON` to the executable path if needed. It regenerates numerical results under `results/` and creates a local `figures/` directory with transcriptomic and metabolic signaling figures. `figures/` is ignored by Git and is not part of the release. Differences in package versions or graphics environments can affect rendering.

## Finding the results

Complete DESeq2 tables, including nonsignificant genes, are in `results/differential_expression/<cell_line>/all_genes.csv`. Normalized counts are in `results/normalized_counts/<cell_line>/normalized_counts.csv`. Complete Hallmark and GO BP enrichment tables, including nonsignificant pathways, and SREBP1-associated module results are in `results/pathway_enrichment/<cell_line>/`. PCA coordinates, explained variance, selected-gene values, and other figure source tables are in `results/plot_data/`. Replace `<cell_line>` with `OVCAR5`, `OVCAR8`, or `ID8`.

## Citation and license

The associated IJMS manuscript is in revision. Its final citation, DOI, and author list should be added after confirmation; no publication DOI is asserted here. For the software repository, use the provisional information in [`CITATION.cff`](CITATION.cff), and cite GEO Series GSE349387 for the sequencing dataset when the record becomes public.

Original analysis code is released under the [MIT License](LICENSE). Third-party gene sets and other source materials retain their own licensing and attribution conditions, described in [`data/gene_sets/README.md`](data/gene_sets/README.md).

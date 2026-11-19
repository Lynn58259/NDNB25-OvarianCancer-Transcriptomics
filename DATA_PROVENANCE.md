# Data and method provenance

The input `data/counts/provider_expression_matrix.tsv` is the unchanged processed expression matrix from the Plasmidsaurus W6XWZQ delivery. It retains provider gene IDs and annotations, fractional counts, CPM, and original sample IDs. `data/metadata/samples.csv` maps those IDs to biological samples, treatments, batch labels, and GEO library names. The RNA-seq analysis uses species-specific count columns; preprocessing rounds them to integers and applies the preserved gene filter before DESeq2. Raw reads are associated with GEO accession GSE349387 and are not included here.

Frozen Hallmark and GO Biological Process gene sets are stored in `data/gene_sets/`; their source release, transformations, attribution, and redistribution terms are documented in [`data/gene_sets/README.md`](data/gene_sets/README.md). `SREBF1_module.rds` and the selected-gene mapping files preserve the prespecified lipogenic analysis inputs. This module is an association-based gene set, not a direct SREBP1-binding measurement.

The release retains the validated count preprocessing, DESeq2 model and contrast, PCA calculation, enrichment ranking and seeds, gene mappings, and plotting code. The batch labels in the analysis metadata should be checked against original laboratory records before publication.

## Release audit provenance (November 18, 2026)

The original provider matrix, all metadata, frozen gene sets, and reference results were preserved byte-for-byte during the audit. SHA-256 hashes of the distributed input files are in `environment/input_sha256.csv`.

The original RNA-seq submission sample list confirms sample numbers 1–18 and their cell-line/treatment/replicate labels. The GEO submission workbook's Metadata sheet, rows 39–56, agrees with all 18 library names and batch numbers. These submission documents do not independently establish passage matching or collection dates. An earlier derived audit records user confirmation of pairing but contains stale `not_confirmed` fields for some pairs. Original laboratory/passage records still need review; the statistical design and batch assignments have not been changed.

The runtime consumes only this repository's processed matrix, sample metadata, frozen RDS gene sets, and selected-gene mapping. External sample records were inspected solely for provenance validation and are not runtime dependencies. No raw read alignment or UMI counting was performed. The distributed provider matrix is SHA-256 identical to the local original delivery; its hash is recorded in `environment/input_sha256.csv`.

GEO public accessibility was not established during this audit. The public GEO URL could not be retrieved by the web tool; this does not prove that the accession is private. No reviewer token is included.

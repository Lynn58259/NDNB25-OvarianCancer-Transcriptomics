# Data and method provenance

The input `data/counts/provider_expression_matrix.tsv` is the unchanged processed expression matrix from the Plasmidsaurus W6XWZQ delivery. It retains provider gene IDs and annotations, fractional counts, CPM, and original sample IDs. `data/metadata/samples.csv` maps those IDs to biological samples, treatments, batch labels, and GEO library names. The RNA-seq analysis uses species-specific count columns; preprocessing rounds them to integers and applies the preserved gene filter before DESeq2. Raw reads are associated with GEO accession GSE349387 and are not included here.

Frozen Hallmark and GO Biological Process gene sets are stored in `data/gene_sets/`; their source release, transformations, attribution, and redistribution terms are documented in [`data/gene_sets/README.md`](data/gene_sets/README.md). `SREBF1_module.rds` and the selected-gene mapping files preserve the prespecified lipogenic analysis inputs. This module is an association-based gene set, not a direct SREBP1-binding measurement.

The release retains the validated count preprocessing, DESeq2 model and contrast, PCA calculation, enrichment ranking and seeds, gene mappings, and plotting code. The batch labels in the analysis metadata should be checked against original laboratory records before publication.

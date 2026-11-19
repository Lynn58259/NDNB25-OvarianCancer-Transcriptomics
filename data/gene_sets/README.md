# Frozen gene sets and ortholog provenance

Source: Molecular Signatures Database (MSigDB), human-source release 2026.1.Hs, retrieved through msigdbr 26.1.0 with babelgene 22.9. Human Hallmark and GO:BP, plus human-source projections to mouse, are stored as named R lists in the four family RDS files. No network call occurs during analysis.

MSigDB copyright (c) 2004-2025 Broad Institute, Inc., Massachusetts Institute of Technology, and Regents of the University of California. These Hallmark and GO collections are provided under Creative Commons Attribution 4.0 International: https://creativecommons.org/licenses/by/4.0/ . See the controlling terms at https://gsea-msigdb.org/gsea/msigdb_license_terms.jsp and resource at https://www.gsea-msigdb.org/gsea/msigdb/ . The repository's MIT license applies to original code and does not replace these third-party terms.

Snapshot transformations: symbols are uppercased for the original matching convention, duplicate members are removed within each set, and gene-set names and membership order are retained. This reproduces the original scripts' processing. Case-preserved mouse symbol evidence is provided separately.

Important reproducibility finding: msigdbr 26.1.0 caches orthologs by species within an R session, using the genes in the first requested collection. Requesting mouse Hallmark before mouse GO in one session truncates the GO mapping. The final mouse GO snapshot was therefore created in a fresh R session, as in the original standalone GO script. All archived selected mouse GO memberships were verified against that snapshot. The normal workflow reads frozen sets and cannot encounter this cache-order issue.

SREBF1_module.rds contains the fixed human genes FASN, ACACA, ACLY, SCD, ELOVL6, FADS1 and FADS2. The ID8 list contains eight projected mouse members. SREBF1_mouse_ortholog_evidence.csv explicitly records SCD -> Scd1 and Scd2 and supporting orthology sources. The fixed module is not a direct SREBF1-binding assay.

selected_gene_mapping.csv freezes the manuscript gene selections and annotation keys, including intentional unavailable human Scd2 rows. It contains no fitted expression values. Gene effects and VST values are regenerated from count-based model outputs.


## Release audit checks

The audit loaded all frozen RDS files and verified their structure, unique pathway names, unique/nonmissing character members, and collection prefixes. Human/mouse Hallmark each contain 50 pathways; human GO BP contains 7,538 and mouse GO BP 7,535. All 10 archived selected mouse GO pathway memberships match the frozen mouse snapshot. The snapshot freezes the historical ortholog projection; it is not a current orthology reannotation.

The [MSigDB licensing page](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) was checked during the audit. The distributed collections are Hallmark and GO BP, with no KEGG or BioCarta collections. Preserve the attribution and transformation notes above when redistributing these snapshots. The processed provider matrix is study data, rather than MIT-licensed original code; the code owner should confirm any provider delivery redistribution conditions.

# Single-cell-transcriptomics-neuroblastoma

The present work gravitates around the implementation of an extensible pipeline of single-cell transcriptomics leveraging RNA-seq counts derived from two MYCN-amplified neuroblastoma cell lines, called KELLY and SK-N-BE(2)-C. The scope of the study is to determine the transcriptome-wide heterogeneity associating with the two different mutational backgrounds harboured by distinct model tumoral cells and evaluate the effectiveness and reuse potential of our generated datsets. 

The flowchart articulates in multiple steps:
1. Knowledge of the datasets and feature-specific variance estimation,
2. Genome-wide mapping of the mean TPM normalized counts,
3. Differential gene expression analysis and differential gene set enrichment analysis,
4. Heterogeneity studies of subpopulations of SK-N-BE(2)-C cells,
5. Estimation and mitigation of the effects of the cell cycle on the exploratory analysis, 
6. Comparison of mean TPMs counts from single-cell RNA-seq and bulk-sequenced transcriptomics data from an expression matrix drawn from the literature.

Single-cell transcriptomics were carried out with the aid of the Seurat package, whereas differential expression analysis was conducted using fgsea combined with gene sets drawn from the molecular signature databases. 






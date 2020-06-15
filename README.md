# Single-cell-transcriptomics-neuroblastoma

The present work gravitates around the implementation of an extensible pipeline of single-cell transcriptomics leveraging RNA-seq counts derived from two MYCN-amplified neuroblastoma cell lines, called KELLY and SK-N-BE(2)-C. The scope of the study is to determine the transcriptome-wide heterogeneity associating with the two different mutational backgrounds harboured by distinct model tumoral cells and evaluate the effectiveness and reuse potential of our generated datsets. 

The flowchart articulates in multiple steps:
1. Knowledge of the datasets and feature-specific variance estimation,
2. Genome-wide mapping of the mean TPM normalized counts,
3. Differential gene expression analysis and differential gene set enrichment analysis,
4. Heterogeneity studies of subpopulations of SK-N-BE(2)-C cells,
5. Estimation and mitigation of the effects of the cell cycle on the exploratory analysis, 
6. Comparison of mean TPMs counts from single-cell RNA-seq and bulk-sequenced transcriptomics data from an expression matrix drawn the 2017 article "Transcriptomic profiling of 39 commonly-used neuroblastoma cell lines" by Harenza et al.

Single-cell transcriptomics were carried out with the aid of the Seurat package, whereas differential expression analysis was conducted using fgsea combined with gene sets drawn from the molecular signature databases. 

##Repository
The repository is organized in three distinct folders:
1. The folder *Data* contains the datasets we utilized during the analysis,
2. The folder *R* stores the R scripts produced for the conduction of the study,
3. The folder *Figures* reports the charts and plots generated by the R scripts during all stages of the analysis,
4. The folder *Tables* reports the tables generated by the R scripts during all stages of the analysis.




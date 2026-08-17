# Description

This repository contains the source code and analysis files for bulk RNA-seq analysis of neuromyelitis optica spectrum disorder (NMOSD) and healthy control samples. The workflow follows bulk RNA-seq processing through transcript quantification, differential expression analysis, and Type I interferon-related gene analysis.

## Dataset

The datasets used in this study involve privacy-sensitive information and are therefore not publicly available.

## Software Versions

The repository uses R scripts and python notebooks. Specific software, package, operating-system, editor, and server versions are not recorded in the repository.

## Analysis Workflow

```mermaid
flowchart TD
    A(["<div style='width:300px;height:90px;text-align:center;'>10 NMOSD bulk RNA-seq<br/>(FASTQ)</div>"])
    B(["<div style='width:300px;height:90px;text-align:center;'>11 healthy control bulk RNA-seq<br/>(FASTQ)</div>"])
    C(["FASTQ quality control<br/>(FastQC / MultiQC)"])
    D(["Adapter trimming<br/>(Cutadapt / Trimmomatic)"])
    E(["Mapping to the reference genome<br/>(STAR / GRCh38)"])
    F(["Transcript identification and gene expression<br/>(RSEM)"])
    G(["Independent cohort<br/>5 NMOSD and 18 healthy controls"])
    H(["Differential expression analysis<br/>(Limma-voom / edgeR / DESeq2)"])
    I(["Type I interferon-related gene set<br/>(MSigDB)"])
    J(["Significant Type I interferon-related genes"])

    A --> C
    B --> C
    C --> D --> E --> F --> H --> I --> J
    G -.-> H

    classDef normal fill:#ffffff,stroke:#555555,stroke-width:1px,color:#111111,font-weight:normal;
    class A,B,C,D,E,F,G,H,I,J normal;
```

The FASTQ quality-control, adapter-trimming, genome-mapping, and RSEM steps are represented in the workflow for completeness. The repository mainly contains the processed RSEM and STAR results used by the downstream analyses.

## Analysis Pipeline and Source Code

| Step | Description | Corresponding code |
|---|---|---|
| Bulk RNA-seq input | The main cohort contains 10 NMOSD and 11 healthy control samples. An independent cohort contains 5 NMOSD and 18 healthy controls. |  |
| FASTQ quality control | Quality assessment of the FASTQ files. |  |
| Adapter trimming | Removal of sequencing adapters. |  |
| Mapping to the reference genome | Alignment to GRCh38. |  |
| Transcript identification and gene expression | Use RSEM and STAR results to construct expression matrices and calculate FPKM/TPM values. | `python/FPKM_TPM.ipynb`, `python/TPM_my_with_external_data.ipynb`, `python/No_3_genes_TPM_my_with_external_data.ipynb`, `python/rna_file.ipynb` |
| Differential expression analysis | Compare NMOSD and healthy control samples in the main and independent cohorts using Limma-voom, edgeR, and DESeq2. | `R/Limma.R`, `R/edgeR.R`, `R/DEseq2.R`, `R/Limma_external.R`, `R/edgeR_external.R`, `R/DEseq2_external.R` |
| Type I interferon-related gene set | Compare DEGs with Type I interferon and hallmark interferon gene sets. | `python/inteferon_geneset.ipynb` |
| Significant Type I interferon-related genes | Integrate DEG results and classify shared up- and down-regulated genes. | `R/DEGs_analysis(3 packages).R`, `python/DEG_up_down.ipynb` |

Additional analysis files include PCA and batch-effect assessment (`R/log2cpm_PCA.R`, `R/TPM_PCA.R`, `R/Batch effect_PCA.R`), hierarchical clustering (`R/hirarchical_clustering.R`), cell-type deconvolution (`R/deconvolution_cell_type_bulk.R`), and functional enrichment analysis (`R/Pathway.R`).

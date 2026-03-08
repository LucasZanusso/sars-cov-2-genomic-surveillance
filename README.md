 # SARS-CoV-2 Genomic Surveillance: Evolutionary Divergence Analysis (2023 vs 2026)

This project implements a complete bioinformatics pipeline for comparative genomic analysis of SARS-

CoV-2, focusing on mutational evolution between 2023 strains and recent 2026 lineages.


# Project Overview

**Data Sourcing**: Automated retrieval of NGS raw data (Fastq) from NCBI SRA, for this project, every sample was extracted from the bioproject PRJNA732685


**Variant Discovery**: Identification of SNPs and Indels relative to the Wuhan reference genome (NC_045512.2).


**Functional Annotation**: Biological impact assessment of mutations using SnpEff.


**Data Visualization**: Comparative genomic tracks and Manhattan plots developed in R.


# Tech Stack

**Bioinformatics**: sra-tools, bwa, samtools, bcftools.


**Annotation**: SnpEff (SARS-CoV-2 database).


**Data Science (R)**: Gviz, ggplot2, VariantAnnotation, Biostrings.


**Workflow**: Bash scripting for pipeline automation and scalability.


# Pipeline Workflow


**Data Acquisition**: Automated download and quality control of SRA accessions.


**Read Mapping**: Alignment against the Wuhan-Hu-1 reference using BWA-MEM.


**Variant Calling**: Germline variant identification with quality filtering (QUAL > 20, DP > 3).


**Genomic Tracks**: Visualization of mutational hotspots across the viral genome.


**Statistical Analysis**: Comparative divergence plots showing mutation accumulation over time.


# Key Results


**Genomic Divergence**


The analysis receals a significant increase in mutational load in 2026 samples compared to 2023 samples.

Hotspots are predominantly located in the Spike (S) Protein region, consistent with positive selection for immune scape


**Genomic Tracks(Gviz)**


<img width="1200" height="1600" alt="comparativo_evolutivo_20231" src="https://github.com/user-attachments/assets/8362688b-a0e6-460b-92d0-b07951e3d1c4" />



<img width="1200" height="1600" alt="comparativo_evolutivo_20261" src="https://github.com/user-attachments/assets/7c2f6f17-9a3f-4045-9ec9-ab76f08c10ba" />


**Variant Quality (Manhattan Plot)**


<img width="3000" height="2400" alt="qualidade_antigas_2023" src="https://github.com/user-attachments/assets/d1132a21-5d36-45bc-bd8b-0117d4e2b47f" />



<img width="3000" height="2400" alt="qualidade_recentes_2026" src="https://github.com/user-attachments/assets/dc4d6587-5639-4a3d-aa9d-536d1bf86708" />


**Total Genomic Divergence**


<img width="3000" height="3000" alt="numero_de_mutações" src="https://github.com/user-attachments/assets/138faa5c-da8b-4b35-9a94-3654e0961449" />


**Evolutionary Insights**

While limited by sample size (n=10), the clear trend in mutation density per kilobase highlights the 

"molecular clock" of SARS-CoV-2, with recent variants exhibiting complex mutational profiles in the 

Receptor Binding Domain (RBD).


 # SARS-CoV-2 Genomic Surveillance: Evolutionary Divergence Analysis (2023 vs 2026)

This project implements a complete bioinformatics pipeline for comparative genomic analysis of SARS-
CoV-2, focusing on mutational evolution between 2023 strains and recent 2026 lineages.


# Project Overview

**Data Sourcing**: Automated retrieval of NGS raw data (Fastq) from NCBI SRA, for this project, every sample was extracted from the bioproject `PRJNA732685`.


**Variant Discovery**: Identification of SNPs and Indels relative to the Wuhan reference genome `(NC_045512.2).`


**Functional Annotation**: Biological impact assessment of mutations using `SnpEff`.


**Data Visualization**: Comparative genomic tracks and Manhattan plots developed in R.


# Tech Stack

**Bioinformatics**: `sra-tools, bwa, samtools, bcftools`.


**Annotation**: `SnpEff` (SARS-CoV-2 database).


**Data Science (R)**:`Gviz, ggplot2, VariantAnnotation, Biostrings`.


**Workflow**: Bash scripting for pipeline automation and scalability.


# Pipeline Workflow


**Data Acquisition**: Automated download and quality control of SRA accessions.


**Read Mapping**: Alignment against the Wuhan-Hu-1 reference using `BWA-MEM`.


**Variant Calling**: Germline variant identification with quality filtering (QUAL > 20, DP > 3).


**Genomic Tracks**: Visualization of mutational hotspots across the viral genome.


**Statistical Analysis**: Comparative divergence plots showing mutation accumulation over time.

# Project Structure & Scripts
Each script in this project has been programmed to be used in sequence, so the output of `00_download_data.sh` will be used as input for `01_qc_batch.sh`, and so on. Therefore, it is necessary to pay attention to the directories created in each step to avoid errors.

### `00_download_data.sh` (SRA Data Retrieval)
This script handles the automated retrieval and organization of raw sequencing data from the **NCBI Sequence Read Archive (SRA)**. It establishes the directory architecture required for the longitudinal comparison.

**Key features:**
* **Automated Workspace Setup:** Initializes the `data/antigo`, `data/recente`, and `logs` directories.
* **Batch Retrieval (`fasterq-dump`):** Efficiently downloads high-throughput sequencing data for 20 accessions.
* **Organized Storage:** Creates individual subdirectories for each `SRR ID`, preventing file collisions.
* **Logging:** Captures output and errors into `logs/download.log`.

**Tools used:** `sra-tools` (`fasterq-dump`).

---

### `01_qc_batch.sh` (Batch QC Pipeline)
This script automates the **Quality Control (QC)** and **Preprocessing** for all sequencing libraries. It ensures that only high-quality reads are used for downstream variant calling.

**Key features:**
* **Directory Management:** Categorizes output in `results/qc/antigo` and `results/qc/recente`.
* **Automated Trimming (`fastp`):** Performs adapter trimming and quality filtering for paired-end reads.
* **Secondary Validation (`FastQC`):** Runs a second round of quality checks on the "cleaned" files.
* **Metadata Aggregation (`MultiQC`):** Summarizes all reports into a single dashboard in `results/qc/multiqc_final`.

**Tools used:** `fastp`, `fastqc`, `multiqc`.

---

### `02_alignment.sh` (Read Mapping & BAM Processing)
This script performs the core task of aligning reads to the **Wuhan-Hu-1 (NC_045512.2)** reference genome using the BWA algorithm.

**Key features:**
* **Reference-Based Alignment (`BWA MEM`):** Executes mapping with multi-threading support (`-t 4`).
* **Read Group Header Addition (`-R`):** Incorporates essential metadata (`ID`, `SM`, `PL:ILLUMINA`) for best practices.
* **BAM Lifecycle Management (`Samtools`):** Handles conversion (SAM to BAM), coordinate sorting, and indexing (`.bai`).
* **Storage Optimization:** Automatically removes intermediate `.sam` and un-sorted `.bam` files to save disk space.

**Tools used:** `bwa`, `samtools`.

---

### `03_variant_calling_.sh` (Genomic Variant Discovery)
This script identifies genetic divergences (SNPs and Indels) between the aligned samples and the reference genome. It implements a robust filtering strategy to ensure the reliability of the identified mutations.

**Key features:**
* **Genomic Likelihood Calculation (`bcftools mpileup`):** Analyzes the sorted `.bam` files to calculate the likelihood of non-reference bases at each genomic position, specifically annotating read depth (`-a DP`).
* **Haploid Variant Calling (`bcftools call`):** Executes the calling model optimized for viral genomes (`--ploidy 1`), focusing only on variant sites (`-v`).
* **Hard Filtering Strategy:** Applies stringent quality thresholds using `bcftools filter` (`QUAL > 20` and `DP > 3`) to eliminate sequencing artifacts and low-confidence calls.
* **Data Optimization:** Generates human-readable `.vcf` files alongside compressed `.vcf.gz` and `tabix` indexes, facilitating high-speed data retrieval for downstream analysis in **R**.
* **Storage Management:** Automatically removes intermediate large `.bcf` files to maintain an efficient project footprint.

**Tools used:** `bcftools`, `htslib` (`bgzip`, `tabix`).

---

### `04_annotation.sh` (Variant Annotation & Impact Assessment)
This script performs the biological interpretation of the identified variants. By integrating the genomic coordinates with the **SARS-CoV-2 (NC_045512.2)** database, it predicts the functional impact of each mutation on the viral proteins.

**Key features:**
* **Functional Prediction (`snpEff`):** Determines the effect of variants on gene expression, identifying synonymous vs. non-synonymous (missense) mutations and stop-gain events.
* **Impact Categorization:** Classifies variants by impact level (HIGH, MODERATE, LOW, or MODIFIER), which is crucial for identifying key changes in the **Spike (S)** and **ORF1ab** regions.
* **Automated Reporting:** Generates comprehensive HTML statistical summaries (`_summary.html`) for each sample, including transitions/transversions ratios and codon change distributions.
* **Resource Optimization:** Utilizes Java heap memory management (`-Xmx4g`) to ensure stable processing of multi-sample datasets.

**Tools used:** `snpEff` (Java-based).

---

---

### `05_genomic_data_analysis.R` (R Data Science & Statistical Assessment)
This comprehensive R script acts as the analytical engine of the project. It integrates genomic variants with biological annotations to perform longitudinal statistical comparisons between 2023 and 2026 SARS-CoV-2 lineages.

**Key Features:**

* **Genomic Track Visualization (`Gviz`):**
    * Constructs high-resolution genomic maps by integrating **TxDb** objects from GFF3 annotations.
    * Generates comparative tracks (2023 vs. 2026) that align mutation hotspots with the viral gene structure (e.g., Spike and ORF1ab proteins).
* **Data Wrangling & Integration (`tidyverse`):**
    * Uses a custom-built function to extract and label variants from multiple VCF files simultaneously.
    * Classifies variants by type (SNP, Insertion, or Deletion) for structural analysis.
* **Statistical Validation:**
    * **Non-Parametric Testing:** Implements the **Wilcoxon Rank Sum Test** to compare mutation counts between periods, accounting for the small sample size ($n=10$).
    * **Categorical Analysis:** Performs **Chi-squared ($\chi^2$) tests** to evaluate if the proportion of mutation types or their location (e.g., inside vs. outside the Spike protein) significantly changed over time.
* **Advanced Data Visualization (`ggplot2`):**
    * **Quality Profiling:** Produces faceted Manhattan-style plots to validate the `QUAL` distribution across all samples.
    * **Evolutionary Divergence:** Generates boxplots and jitter plots illustrating the "Genetic Distance" from the Wuhan reference strain.
    * **Proportional Analysis:** Visualizes the shift in variant types using stacked percentage bar charts, integrated with p-values from statistical tests.

**Libraries used:** `Gviz`, `VariantAnnotation`, `GenomicFeatures`, `tidyverse`, `ggplot2`, `Biostrings`.

---
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

# How to reproduce:
1. Download the scripts.sh\
2. Install dependencies (Conda environment recommended)\
3. Enter the desired SRR identifiers in the download_data.sh script\
4. Run the scripts in sequence

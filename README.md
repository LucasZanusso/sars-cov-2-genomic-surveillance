 # SARS-CoV-2 Genomic Surveillance: Evolutionary Divergence Analysis (2023 vs 2026)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Language: Bash](https://img.shields.io/badge/Language-Bash-green.svg)
![Language: R](https://img.shields.io/badge/Language-R-blue.svg)
![Platform: Linux](https://img.shields.io/badge/Platform-Linux-lightgrey.svg)

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

| Step | Script | Tool(s) | Output |
|------|--------|---------|--------|
| 1. Data Retrieval | `00_download_data.sh` | sra-tools | `.fastq` |
| 2. Quality Control | `01_qc_batch.sh` | fastp, FastQC, MultiQC | trimmed `.fastq`, QC reports |
| 3. Alignment | `02_alignment.sh` | BWA-MEM, Samtools | sorted `.bam` |
| 4. Variant Calling | `03_variant_calling.sh` | BCFtools | `.vcf.gz` |
| 5. Annotation | `04_annotation.sh` | SnpEff | annotated `.vcf`, HTML report |
| 6. Analysis | `05_genomic_data_analysis.R` | R/Bioconductor | figures, statistics |


# Project Structure & Scripts
Each script in this project has been programmed to be used in sequence, so the output of `00_download_data.sh` will be used as input for `01_qc_batch.sh`, and so on. Therefore, it is necessary to pay attention to the directories created in each step to avoid errors.

### Prerequisites
`environment.yml` — Conda environment with all dependencies and versions
`samples.txt` — SRR accessions used in this analysis

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

### `03_variant_calling.sh` (Genomic Variant Discovery)
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


The analysis reveals a significant increase in mutational load in 2026 samples compared to 2023 samples.

Hotspots are predominantly located in the Spike (S) Protein region, consistent with positive selection for immune escape


**Genomic Tracks(Gviz)**

| 2023 Variants | 2026 Variants |
|:---:|:---:|
| <img width="500" src="https://github.com/user-attachments/assets/89dc947e-fe16-42ed-9628-80d35c44efcb" /> | <img width="500" src="https://github.com/user-attachments/assets/9b80e562-893c-45eb-8945-1a6d4b779d47" />|

**Figure 1.** Genomic tracks showing variant distribution across the SARS-CoV-2 genome 
(NC_045512.2). Each row represents one sample. Colored marks indicate SNPs and Indels 
mapped against the annotated gene structure. The 2026 lineages show a visibly higher 
mutation density, particularly in the Spike (S) gene region (positions ~21,563–25,384).

---

### Variant Quality Distribution (Manhattan Plot)

| 2023 Samples | 2026 Samples |
|:---:|:---:|
| <img width="500" src="https://github.com/user-attachments/assets/d350fe3a-3227-4751-98de-197983f8c629" />| <img width="500" src="https://github.com/user-attachments/assets/6769a018-a2e3-4ff1-84b8-9aaf496cf1b6" /> |

**Figure 2.** QUAL score distribution per genomic position across all samples. 
Only variants passing the hard filter (QUAL > 20, DP > 3) are shown. 
Both cohorts display comparable quality profiles, supporting the validity of 
the downstream comparative analysis.

---

### Total Genomic Divergence

<img width="700" src="https://github.com/user-attachments/assets/c26d32fc-7934-460a-bc5c-719bc7155207" />

**Figure 3.** Total mutation count per sample relative to the Wuhan-Hu-1 reference. 
2026 lineages accumulate significantly more variants than 2023 strains 
(Wilcoxon rank-sum test, p < 0.05), consistent with the ongoing molecular 
evolution of SARS-CoV-2 under immune selection pressure.

---

## Evolutionary Insights

The integrative analysis of mutational load, variant composition, and functional annotation 
reveals consistent molecular signatures of adaptive evolution in 2026 SARS-CoV-2 lineages.

**Directional increase in mutational load.** 2026 samples accumulated significantly more 
variants relative to the Wuhan-Hu-1 reference than 2023 strains (Wilcoxon rank-sum test, 
W = 0, p = 0.0004), with a W statistic of zero indicating complete separation between groups — 
every 2026 sample harbored more mutations than any 2023 sample.

**Qualitative shift in variant composition.** Beyond total count, the proportion of variant 
types differed significantly between periods (χ² = 6.04, df = 1, p = 0.014), indicating that 
the mutational landscape is not only larger but qualitatively distinct in 2026 — consistent 
with directional selection rather than neutral drift.

**Positive selection concentrated in the Spike protein.** The proportion of variants mapping 
to the Spike (S) gene was significantly higher in 2026 samples (χ² = 6.32, df = 1, p = 0.012). 
Notably, 2026 lineages carry ~67 moderate-impact variants per sample in the Spike alone, 
compared to ~2 in 2023 samples — a ~33-fold increase concentrated at the primary antigenic 
target of the virus.

**Fixed immune escape mutations in the RBD.** All 10 (10/10) 2026 samples share a conserved 
set of non-synonymous substitutions in the Receptor Binding Domain, including positions 
previously associated with immune escape and ACE2 binding affinity: S477N, T478K, F486L/S, 
F456L, L455S, and K417N. The co-occurrence of F456L + L455S + F486L/S is a hallmark of 
JN.1-derived lineages (e.g., KP.x, XEC), consistent with the variants circulating globally 
during this period. The P681R substitution at the furin cleavage site — a Delta lineage 
signature associated with increased fusogenicity — was also fixed across all 2026 samples, 
suggesting convergent evolution at this functionally critical position.

> **Note:** These findings are exploratory and limited by sample size (n = 10 per group). 
> Variant interpretation is based on SnpEff functional annotation against NC_045512.2. 
> Lineage assignment would require phylogenetic analysis (e.g., Nextclade/PANGOLIN).

# How to Reproduce

1. Clone this repository:
   `git clone https://github.com/LucasZanusso/sars-cov-2-genomic-surveillance`
2. Create and activate the Conda environment:
   `conda env create -f environment.yml && conda activate sars-cov2-surveillance`
3. Edit `samples.txt` with your desired SRR accessions, or use the provided ones
4. Run the scripts in sequence: `bash 00_download_data.sh`, and so on

> **Note:** The `results/` folder contains a MultiQC report, example VCFs, and SnpEff annotations for reference.

# =============================================================================
# 06_spike_mutation_analysis.R
# Spike Protein Mutation Analysis - SARS-CoV-2 2023 vs 2026
# 
# Outputs:
#   - results/figures/spike_variant_composition.png  (barplot missense vs synonymous)
#   - results/tables/spike_comparative_table.tsv     (comparative table 2023 vs 2026)
#   - results/tables/fixed_mutations_2026.tsv        (mutations fixed in 10/10 samples)
# =============================================================================

library(tidyverse)

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR      <- "/home/lucas/covid_project_final"
ANNOT_DIR     <- file.path(BASE_DIR, "results/annotation")
FIGURES_DIR   <- file.path(BASE_DIR, "results/figures")
TABLES_DIR    <- file.path(BASE_DIR, "results/tables")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)

# Known immune escape positions in the Spike RBD (literature-curated)
ESCAPE_POSITIONS <- c(417, 445, 452, 455, 456, 477, 478, 484, 486, 490, 493, 494, 
                      501, 505, 681)

# =============================================================================
# FUNCTION: Parse _summary.genes.txt files and extract Spike data
# =============================================================================
parse_spike_summary <- function(filepath, sample_id, periodo) {
  
  # Lê sem assumir número de colunas
  df <- read_tsv(filepath, comment = "#", show_col_types = FALSE,
                 col_names = FALSE)
  
  # Filtra só a linha do Spike (transcript principal)
  df_spike <- df %>%
    filter(X1 == "S", X3 == "GU280_gp02")
  
  if (nrow(df_spike) == 0) {
    message("  WARNING: No Spike data found in ", sample_id)
    return(tibble(sample_id = sample_id, periodo = periodo,
                  missense = NA_real_, synonymous = NA_real_,
                  impact_low = NA_real_, impact_mod = NA_real_))
  }
  
  # Pega a última linha (algumas amostras têm duplicatas)
  row <- df_spike %>% slice_tail(n = 1)
  
  # As colunas de impacto e efeito são as últimas 9 colunas sempre
  ncols <- ncol(row)
  
  tibble(
    sample_id  = sample_id,
    periodo    = periodo,
    impact_low = as.numeric(row[[ncols - 8]]),
    impact_mod = as.numeric(row[[ncols - 7]]),
    missense   = as.numeric(row[[ncols - 2]]),
    synonymous = as.numeric(row[[ncols - 1]])
  )
}
# =============================================================================
# FUNCTION: Parse VCF and extract Spike amino acid changes
# =============================================================================

parse_spike_summary <- function(filepath, sample_id, periodo) {
  
  df <- read_tsv(filepath, comment = "#", show_col_types = FALSE,
                 col_names = FALSE)
  
  # Arquivo vazio — retorna NA
  if (nrow(df) == 0 || ncol(df) == 0) {
    message("  WARNING: Empty file for ", sample_id)
    return(tibble(sample_id = sample_id, periodo = periodo,
                  missense = NA_real_, synonymous = NA_real_,
                  impact_low = NA_real_, impact_mod = NA_real_))
  }
  
  df_spike <- df %>%
    filter(X1 == "S", X3 == "GU280_gp02")
  
  if (nrow(df_spike) == 0) {
    message("  WARNING: No Spike data found in ", sample_id)
    return(tibble(sample_id = sample_id, periodo = periodo,
                  missense = NA_real_, synonymous = NA_real_,
                  impact_low = NA_real_, impact_mod = NA_real_))
  }
  
  row <- df_spike %>% slice_tail(n = 1)
  ncols <- ncol(row)
  
  tibble(
    sample_id  = sample_id,
    periodo    = periodo,
    impact_low = as.numeric(row[[ncols - 8]]),
    impact_mod = as.numeric(row[[ncols - 7]]),
    missense   = as.numeric(row[[ncols - 2]]),
    synonymous = as.numeric(row[[ncols - 1]])
  )
}
# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading summary.genes.txt files...")

summary_files <- tibble(
  filepath  = c(
    list.files(file.path(ANNOT_DIR, "antigo"),  pattern = "summary.genes.txt",
               recursive = TRUE, full.names = TRUE),
    list.files(file.path(ANNOT_DIR, "recente"), pattern = "summary.genes.txt",
               recursive = TRUE, full.names = TRUE)
  )
) %>%
  mutate(
    sample_id = basename(dirname(filepath)),
    periodo   = ifelse(str_detect(filepath, "/antigo/"), "2023", "2026")
  )
for (i in seq_len(nrow(summary_files))) {
  df_test <- read_tsv(summary_files$filepath[i], comment = "#", 
                      show_col_types = FALSE, col_names = FALSE)
  cat(i, summary_files$sample_id[i], "| nrow:", nrow(df_test), 
      "| ncol:", ncol(df_test), "\n")
}
df_summary <- pmap_dfr(
  list(summary_files$filepath, summary_files$sample_id, summary_files$periodo),
  parse_spike_summary
)

message("Loading VCF files for amino acid changes...")

vcf_files <- tibble(
  filepath = c(
    list.files(file.path(ANNOT_DIR, "antigo"),  pattern = "\\.vcf$",
               recursive = TRUE, full.names = TRUE),
    list.files(file.path(ANNOT_DIR, "recente"), pattern = "\\.vcf$",
               recursive = TRUE, full.names = TRUE)
  )
) %>%
  mutate(
    sample_id = basename(dirname(filepath)),
    periodo   = ifelse(str_detect(filepath, "/antigo/"), "2023", "2026")
  )

df_vcf <- pmap_dfr(
  list(vcf_files$filepath, vcf_files$sample_id, vcf_files$periodo),
  parse_spike_vcf
)

# =============================================================================
# OUTPUT 1: BARPLOT — Missense vs Synonymous per period
# =============================================================================

message("Generating barplot: missense vs synonymous...")

df_bar <- df_summary %>%
  group_by(periodo) %>%
  summarise(
    Missense   = mean(missense,   na.rm = TRUE),
    Synonymous = mean(synonymous, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Missense, Synonymous),
               names_to  = "variant_type",
               values_to = "mean_count")

p_bar <- ggplot(df_bar, aes(x = periodo, y = mean_count, fill = variant_type)) +
  geom_col(position = "dodge", width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(mean_count, 1)),
            position = position_dodge(width = 0.6),
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(
    values = c("Missense" = "#E63946", "Synonymous" = "#457B9D"),
    name   = "Variant Type"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Spike Protein: Mean Variant Count per Period",
    subtitle = "Missense (non-synonymous) vs Synonymous variants | n = 10 per group",
    x        = "Period",
    y        = "Mean number of variants per sample",
    caption  = "Source: SnpEff annotation against NC_045512.2"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position = "top",
    axis.line     = element_line(color = "grey30"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

ggsave(file.path(FIGURES_DIR, "spike_variant_composition.png"),
       plot = p_bar, width = 7, height = 5, dpi = 300)
message("  Saved: results/figures/spike_variant_composition.png")

# =============================================================================
# OUTPUT 2: COMPARATIVE TABLE — 2023 vs 2026 per sample
# =============================================================================

message("Generating comparative table...")

df_table <- df_summary %>%
  select(sample_id, periodo, missense, synonymous, impact_low, impact_mod) %>%
  arrange(periodo, sample_id)

# Summary rows
df_summary_row <- df_table %>%
  group_by(periodo) %>%
  summarise(
    sample_id  = paste0("MEAN (", first(periodo), ")"),
    missense   = round(mean(missense,   na.rm = TRUE), 1),
    synonymous = round(mean(synonymous, na.rm = TRUE), 1),
    impact_low = round(mean(impact_low, na.rm = TRUE), 1),
    impact_mod = round(mean(impact_mod, na.rm = TRUE), 1),
    .groups = "drop"
  )

df_table_final <- bind_rows(df_table, df_summary_row) %>%
  rename(
    Sample        = sample_id,
    Period        = periodo,
    Missense      = missense,
    Synonymous    = synonymous,
    Impact_LOW    = impact_low,
    Impact_MODERATE = impact_mod
  )

write_tsv(df_table_final, file.path(TABLES_DIR, "spike_comparative_table.tsv"))
message("  Saved: results/tables/spike_comparative_table.tsv")

# =============================================================================
# OUTPUT 3: FIXED MUTATIONS — present in 10/10 2026 samples
# =============================================================================

message("Identifying fixed mutations in 2026 (10/10 samples)...")

n_recente <- df_vcf %>% filter(periodo == "2026") %>% pull(sample_id) %>% n_distinct()

df_fixed <- df_vcf %>%
  filter(periodo == "2026", effect == "missense_variant") %>%
  group_by(aa_change) %>%
  summarise(
    n_samples   = n_distinct(sample_id),
    frequency   = round(n_samples / n_recente, 2),
    .groups     = "drop"
  ) %>%
  filter(n_samples == n_recente) %>%
  mutate(
    position = as.integer(str_extract(aa_change, "[0-9]+")),
    immune_escape = position %in% ESCAPE_POSITIONS
  ) %>%
  arrange(position) %>%
  rename(
    AA_Change     = aa_change,
    N_Samples     = n_samples,
    Frequency     = frequency,
    Position      = position,
    Immune_Escape = immune_escape
  )

write_tsv(df_fixed, file.path(TABLES_DIR, "fixed_mutations_2026.tsv"))

# Print summary to console
cat("\n========================================\n")
cat("FIXED MUTATIONS IN 2026 (", n_recente, "/", n_recente, " samples)\n")
cat("========================================\n")
print(df_fixed, n = Inf)

n_escape <- sum(df_fixed$Immune_Escape)
cat("\nTotal fixed missense mutations:", nrow(df_fixed))
cat("\nAt known immune escape positions:", n_escape, 
    paste0("(", round(n_escape/nrow(df_fixed)*100, 1), "%)"), "\n")
cat("========================================\n\n")

message("Done. All outputs saved to results/figures/ and results/tables/")

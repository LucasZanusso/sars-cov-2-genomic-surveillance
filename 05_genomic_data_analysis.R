#==============================================================================#
# 1. CONFIGURAÇÕES E DEPENDÊNCIAS
#==============================================================================#

# Bibliotecas
library(tidyverse)
library(VariantAnnotation)
library(GenomicFeatures)
library(Gviz)
library(txdbmaker)
library(scales)

# Configurações globais
base_path <- "/home/lucas/covid_project_final"
setwd(base_path)
options(ucscChromosomeNames = FALSE)

# Definição de Cores e ROI
cores_projeto <- c("2023" = "#999999", "2026" = "#002060")
pos_spike <- c(start = 21563, end = 25384)

#==============================================================================#
# 2. TABELA MESTRE DE AMOSTRAS (Metadados)
#==============================================================================#

# Criamos um dataframe que centraliza tudo o que sabemos sobre as amostras
df_amostras <- tibble(
  id = c("SRR24510773", "SRR24510774", "SRR24510775", "SRR24510776", "SRR24510777", 
  "SRR24510778", "SRR24510780", "SRR24510779", "SRR24510782",
         "SRR37461704", "SRR37461705", "SRR37461706", "SRR37461707", "SRR37461708",
  "SRR37461709", "SRR37461710", "SRR37461711", "SRR37461712"),
  periodo = rep(c("2023", "2026"), each = 9),
  pasta   = rep(c("antigo", "recente"), each = 9)
) %>%
  mutate(
    vcf_path = glue::glue("{base_path}/results/vcf/{pasta}/{id}/{id}_filtered.vcf")
  )

#==============================================================================#
# 3. INFRAESTRUTURA GENÔMICA
#==============================================================================#
# Definindo a ROI (Region of Interest)
# Aqui você pode definir o genoma inteiro ou uma região específica
# Exemplo: Genoma completo do SARS-CoV-2 (aprox. 29.9kb)
roi <- GRanges(
  seqnames = "NC_045512.2", 
  ranges = IRanges(start = 1, end = 29903))
  
  
gff_file <- file.path(base_path, "ref/NC_045512.2.gff3")
txdb <- makeTxDbFromGFF(gff_file, format = "gff3")
seqlevels(txdb, pruning.mode="coarse") <- "NC_045512.2"

# Track de referência (Wuhan)
geneTrack <- GeneRegionTrack(
  txdb, 
  genome = "covid", 
  chromosome = "NC_045512.2", 
  name = "Genes (Wuhan)",
  fill = "orange",
  transcriptAnnotation = "symbol",
  background.title = "brown"
)

genomeAxisTrack <- GenomeAxisTrack()


#==============================================================================#
# 4. CRIAÇÃO AUTOMATIZADA DE TRACKS (Gviz)
#==============================================================================#

# Definimos uma função para criar cada track individualmente
# Isso facilita ajustes estéticos futuros em um único lugar
criar_track_vcf <- function(vcf_path, id, periodo, roi_obj) {
  
  # Leitura e filtro
  vcf_temp <- readVcf(vcf_path, genome = "covid")
  vcf_roi_temp <- subsetByOverlaps(vcf_temp, roi_obj)
  
  # Seleção da cor baseada no período definido na nossa tabela mestre
  cor_fill <- cores_projeto[periodo]
  
  AnnotationTrack(
    start = start(rowRanges(vcf_roi_temp)),
    end   = end(rowRanges(vcf_roi_temp)),
    chromosome = "NC_045512.2",
    name = id,
    shape = "box",
    fill = cor_fill,
    col = "transparent"
  )
}

# Criamos todos os tracks de uma vez e guardamos na própria tabela df_amostras
# Nota: certifique-se que o objeto 'roi' já está definido no seu ambiente
df_amostras <- df_amostras %>%
  mutate(
    track = pmap(list(vcf_path, id, periodo), ~criar_track_vcf(..1, ..2, ..3, roi))
  )

#==============================================================================#
# 5. GERAÇÃO DOS GRÁFICOS COMPARATIVOS
#==============================================================================#

# Função auxiliar para plotar por grupo, evitando repetir o código do plotTracks
gerar_plot_comparativo <- function(dados_grupo, titulo, cor_fundo, arquivo_saida) {
  png(arquivo_saida, width = 1200, height = 1600, res = 150)
  
  # Monta a lista de tracks: Eixo + Genes + Tracks das Amostras
  lista_tracks <- c(list(genomeAxisTrack, geneTrack), dados_grupo$track)
  
  plotTracks(
    lista_tracks,
    from = 21000, to = 29903,
    chromosome = "NC_045512.2",
    main = titulo,
    cex.main = 1,
    background.title = cor_fundo
  )
  dev.off()
}

# Gerar os dois PDFs/PNGs filtrando a tabela mestre
df_amostras %>%
  filter(periodo == "2023") %>%
  gerar_plot_comparativo("Amostras Antigas (2023)", "grey30", "comparativo_2023.png")

df_amostras %>%
  filter(periodo == "2026") %>%
  gerar_plot_comparativo("Amostras Recentes (2026)", "#002060", "comparativo_2026.png")


#==============================================================================#
# 6. EXTRAÇÃO UNIFICADA DE DADOS (Tidy Data)
#==============================================================================#

extrair_dados_vcf <- function(path, id, periodo) {
  vcf_tmp <- readVcf(path, genome = "covid")
  vr_tmp <- rowRanges(vcf_tmp)
  
  # Lógica de tipos de variante otimizada
  ref_seq <- as.character(ref(vcf_tmp))
  alt_seq <- sapply(alt(vcf_tmp), function(x) as.character(x[1]))
  
  tipo_var <- case_when(
    nchar(ref_seq) == nchar(alt_seq) ~ "SNP",
    nchar(ref_seq) < nchar(alt_seq)  ~ "INS",
    TRUE                            ~ "DEL"
  )
  
  tibble(
    pos = start(vr_tmp),
    qual = qual(vcf_tmp),
    amostra = id,
    tipo = tipo_var,
    periodo = periodo
  ) %>%
    # Já adicionamos a classificação da Spike aqui na extração
    mutate(localizacao = if_else(pos >= pos_spike["start"] & pos <= pos_spike["end"], 
                                 "Spike", "Outros"))
}

# Aplicamos a função a todas as amostras da nossa tabela mestre
df_total <- pmap_dfr(list(df_amostras$vcf_path, df_amostras$id, df_amostras$periodo), 
                     extrair_dados_vcf)

#==============================================================================#
# 7. VISUALIZAÇÃO DE QUALIDADE (ggplot2)
#==============================================================================#

# Criamos uma função para evitar repetir o código do ggplot para 2023 e 2026
plotar_qualidade <- function(df, ano, cor_ponto) {
  df %>%
    filter(periodo == ano) %>%
    ggplot(aes(x = pos, y = qual)) +
    geom_point(alpha = 0.5, color = cor_ponto) +
    facet_wrap(~amostra, ncol = 1) + 
    labs(title = paste("Qualidade das Variantes:", ano),
         subtitle = "Referência: Wuhan (NC_045512.2)",
         x = "Posição Genómica", y = "Qualidade (QUAL)") +
    theme_minimal()
}

# Gerar e guardar
p_2023 <- plotar_qualidade(df_total, "2023", cores_projeto["2023"])
p_2026 <- plotar_qualidade(df_total, "2026", cores_projeto["2026"])

ggsave("qualidade_2023.png", p_2023, width = 10, height = 8)
ggsave("qualidade_2026.png", p_2026, width = 10, height = 8)

#==============================================================================#
# 8. ANÁLISE ESTATÍSTICA E DIVERGÊNCIA
#==============================================================================#

# Resumo por amostra
resumo_estatistico <- df_total %>%
  group_by(amostra, periodo) %>%
  summarise(
    Total_Mutacoes = n(),
    Qualidade_Media = mean(qual, na.rm = TRUE),
    SNPs = sum(tipo == "SNP"),
    Indels = sum(tipo %in% c("INS", "DEL")),
    Mut_por_kb = n() / 29.9,
    .groups = "drop"
  )

# Teste de Wilcoxon para o número de mutações entre períodos
teste_mutacoes <- wilcox.test(Total_Mutacoes ~ periodo, data = resumo_estatistico)

# Tabela de Contingência (Tipo de Variante vs Período)
tabela_tipo <- table(df_total$periodo, df_total$tipo)
teste_qui_tipo <- chisq.test(tabela_tipo)

# Tabela de Contingência (Localização Spike vs Período)
tabela_spike <- table(df_total$periodo, df_total$localizacao)
teste_qui_spike <- chisq.test(tabela_spike)

#==============================================================================#
# 9. PLOT DE DIVERGÊNCIA FINAL
#==============================================================================#
ggplot(resumo_estatistico, aes(x = periodo, y = Total_Mutacoes, fill = periodo)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) 
  geom_jitter(width = 0.1, size = 3, aes(color = periodo)) +
  scale_fill_manual(values = cores_projeto) +
  scale_color_manual(values = cores_projeto) +
  labs(title = "Divergência Genética: 2023 vs 2026",
       subtitle = paste("p-valor (Wilcoxon):", round(teste_mutacoes$p.value, 4)),
       x = "Período", y = "Nº de Mutações (Distância de Wuhan)") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("divergencia_final.png", width = 8, height = 6)





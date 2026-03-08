


#==================================Parte_1=====================================#

setwd('/home/lucas/teste')

library(VariantAnnotation)
library(GenomicFeatures)
library(GenomicRanges)
library(Gviz)
library(Biostrings)
library(txdbmaker) 
library(VariantAnnotation)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GenomeInfoDb)
library(ggplot2)
library(tidyverse)

amostras <- c("SRR24510773" = "/home/lucas/teste/results/vcf/antigo/SRR24510773/SRR24510773_filtered.vcf",
              "SRR24510774" = "/home/lucas/teste/results/vcf/antigo/SRR24510774/SRR24510774_filtered.vcf",
              "SRR24510775" = "/home/lucas/teste/results/vcf/antigo/SRR24510775/SRR24510775_filtered.vcf",
              "SRR24510776" = "/home/lucas/teste/results/vcf/antigo/SRR24510776/SRR24510776_filtered.vcf",
              "SRR24510777" = "/home/lucas/teste/results/vcf/antigo/SRR24510777/SRR24510777_filtered.vcf",
              "SRR37461704"   = "/home/lucas/teste/results/vcf/recente/SRR37461704/SRR37461704_filtered.vcf",
              "SRR37461705"   = "/home/lucas/teste/results/vcf/recente/SRR37461705/SRR37461705_filtered.vcf",
              "SRR37461706"   = "/home/lucas/teste/results/vcf/recente/SRR37461706/SRR37461706_filtered.vcf",
              "SRR37461707"   = "/home/lucas/teste/results/vcf/recente/SRR37461707/SRR37461707_filtered.vcf",
              "SRR37461708"   = "/home/lucas/teste/results/vcf/recente/SRR37461708/SRR37461708_filtered.vcf")

options(ucscChromosomeNames = FALSE)
snp_tracks_list <- list()

gff_file <- "/home/lucas/teste/ref/NC_045512.2.gff3"
txdb <- makeTxDbFromGFF(gff_file, format = "gff3")


seqlevels(txdb, pruning.mode="coarse") <- "NC_045512.2"

#Criar o GeneRegionTrack (o geneTrack)
geneTrack <- GeneRegionTrack(
  txdb, 
  genome = "covid", 
  chromosome = "NC_045512.2", 
  name = "Genes (Wuhan)",
  fill = "orange",
  transcriptAnnotation = "symbol",
  background.title = "brown"
)


for(i in 1:length(amostras)) {
  nome_id <- names(amostras)[i]
  vcf_path <- amostras[i]
  
  #Leitura e subset por ROI
  vcf_temp <- readVcf(vcf_path, genome = "covid")
  vcf_roi_temp <- subsetByOverlaps(vcf_temp, roi)
  
  #Cores: Cinza para 2023, Azul Escuro para 2026
  cor_track <- ifelse(grepl("Antiga", nome_id), "#999999", "#002060")
  
  snp_tracks_list[[nome_id]] <- AnnotationTrack(
    start = start(rowRanges(vcf_roi_temp)),
    end   = end(rowRanges(vcf_roi_temp)),
    chromosome = "NC_045512.2",
    name = nome_id,
    shape = "box",
    fill = cor_track,
    col = "transparent" # Deixa o visual mais limpo
  )
}

# Criar lista apenas com as Antigas
tracks_antigas <- snp_tracks_list[grepl("SRR2451077", names(snp_tracks_list))]

# Criar lista apenas com as Recentes
tracks_recentes <- snp_tracks_list[grepl("SRR3746170", names(snp_tracks_list))]

# Configura a tela para ter 2 linhas e 1 coluna
par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))

# Gráfico 1: Antigas (Cinza)
png("comparativo_evolutivo_2023.png", width = 1200, height = 1600, res = 150)
plotTracks(
  c(list(genomeAxisTrack, geneTrack), tracks_antigas),
  from = 21000, to = 29903,
  chromosome = "NC_045512.2",
  main = "Cargas Mutacionais - Amostras Antigas (2023)",
  cex.main = 1,
  background.title = "grey30"
)

# Gráfico 2: Recentes (Azul)
png("comparativo_evolutivo_2026.png", width = 1200, height = 1600, res = 150)
plotTracks(
  c(list(genomeAxisTrack, geneTrack), tracks_recentes),
  from = 21000, to = 29903,
  chromosome = "NC_045512.2",
  main = "Cargas Mutacionais - Amostras Recentes (2026)",
  cex.main = 1,
  background.title = "#002060"
)


#==================================Parte_2=====================================#


# Função para extrair dados e rotular
extrair_dados_vcf <- function(path, label) {
  vcf_tmp <- readVcf(path, genome = "covid")
  vr_tmp <- rowRanges(vcf_tmp)
  
  # Determina tipo de variante
  ref_seq <- as.character(ref(vcf_tmp))
  alt_seq <- sapply(alt(vcf_tmp), function(x) as.character(x[1]))
  type <- ifelse(nchar(ref_seq) == nchar(alt_seq), "SNP",
                 ifelse(nchar(ref_seq) < nchar(alt_seq), "INS", "DEL"))
  
  data.frame(
    pos = start(vr_tmp),
    qual = qual(vcf_tmp),
    amostra = label,
    tipo = type,
    periodo = ifelse(grepl("SRR2451077", label), "2023", "2026")
  )
}

# Une os dados das amostras
df_total <- map2_df(amostras, names(amostras), extrair_dados_vcf)

# Filtrar apenas o grupo antigo
df_antigo <- df_total %>% filter(periodo == "2023")

# Plot para as Antigas
ggplot(df_antigo, aes(x = pos, y = qual)) +
  geom_point(alpha = 0.5, color = "gray") +
  facet_wrap(~amostra, ncol = 1) + 
  labs(title = "Distribuição de Qualidade: Amostras Antigas (2023)",
       subtitle = "Variantes identificadas em relação à referência de Wuhan",
       x = "Posição Genômica (NC_045512.2)", 
       y = "Qualidade (QUAL)") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

# Para salvar:
ggsave("qualidade_antigas_2023.png", width = 10, height = 8)

# Filtrar apenas o grupo recente
df_recente <- df_total %>% filter(periodo == "2026")

# Plot para as Recentes
ggplot(df_recente, aes(x = pos, y = qual)) +
  geom_point(alpha = 0.5, color = "#002060") +
  facet_wrap(~amostra, ncol = 1) + 
  labs(title = "Distribuição de Qualidade: Amostras Recentes (2026)",
       subtitle = "Variantes identificadas em relação à referência de Wuhan",
       x = "Posição Genômica (NC_045512.2)", 
       y = "Qualidade (QUAL)") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

# Para salvar:
ggsave("qualidade_recentes_2026.png", width = 10, height = 8)

#==========================Divergência_genética================================#

# 1. Calcular o número de mutações por amostra
distancia_dados <- df_total %>%
  group_by(amostra, periodo) %>%
  summarise(distancia = n()) 

# 2. Teste de Correlação ou ANOVA simples
# Como temos apenas dois grupos, vamos usar um Boxplot com os pontos individuais
ggplot(distancia_dados, aes(x = periodo, y = distancia, fill = periodo)) +
  geom_boxplot(alpha = 0.3) +
  geom_jitter(width = 0.1, size = 3) +
  labs(title = "Divergência Genética Total (Wuhan vs. Amostras)",
       subtitle = "Acúmulo de mutações ao longo do tempo (Diferença Biológica)",
       x = "Período",
       y = "Número de Mutações (Distância da Referência)") +
  theme_minimal()
ggsave('numero_de_mutações.png', width = 10, height = 10)


#=======================Testes_estatísticos====================================#


contagem <- df_total %>% group_by(amostra, periodo) %>% summarise(n = n())

# Teste t (se os dados forem normais) ou Wilcoxon (mais seguro para n pequeno)
wilcox.test(n ~ periodo, data = contagem)


resumo_estatistico <- df_total %>%
  group_by(amostra, periodo) %>%
  summarise(
    Total_Mutacoes = n(),
    Qualidade_Media = mean(qual),
    SNPs = sum(tipo == "SNP"),
    Indels = sum(tipo %in% c("INS", "DEL")),
    Mut_por_kb = n() / 29.9 # Tamanho do genoma em kb
  )

print(resumo_estatistico)


# Exemplo de tabela de contingência
tabela <- matrix(c(40, 5,  # 2020: 40 SNPs, 5 Indels
                   120, 35), # 2024: 120 SNPs, 35 Indels
                 nrow = 2, byrow = TRUE)
chisq.test(tabela)

#================================Contingencia==================================#

# 1. Agrupar os dados por período e tipo de variante
tabela_contingencia <- df_total %>%
  mutate(categoria_tipo = ifelse(tipo == "SNP", "SNP", "INDEL")) %>%
  group_by(periodo, categoria_tipo) %>%
  summarise(contagem = n()) %>%
  spread(categoria_tipo, contagem, fill = 0) %>%
  column_to_rownames("periodo")

# Visualizar a tabela
print(tabela_contingencia)

# 2. Realizar o Teste de Qui-quadrado
teste_qui2 <- chisq.test(tabela_contingencia)

# 3. Ver o resultado
print(teste_qui2)  

ggplot(df_total, aes(x = periodo, fill = tipo)) +
  geom_bar(position = "fill") + # "fill" mostra a proporção (0 a 100%)
  labs(title = "Proporção de Tipos de Variantes por Período",
       subtitle = paste("p-valor (Qui-quadrado):", round(teste_qui2$p.value, 4)),
       x = "Ano da Amostra",
       y = "Proporção",
       fill = "Tipo") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

# Crie uma coluna para ver se a mutação está na Spike (posições ~21563 a 25384)
df_total <- df_total %>%
  mutate(is_spike = ifelse(pos >= 21563 & pos <= 25384, "Spike", "Outros"))

# Tabela de contingência: Ano vs Localização
tabela_spike <- table(df_total$periodo, df_total$is_spike)
chisq.test(tabela_spike)



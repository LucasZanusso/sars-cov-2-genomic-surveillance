#!/bin/bash

# Criar pasta geral de resultados de QC
mkdir -p results/qc/antigo
mkdir -p results/qc/recente

# Função para processar uma categoria (antigo ou recente)
process_category() {
    local categoria=$1  # 'antigo' ou 'recente'
    
    # Percorre cada subpasta de ID dentro de data/categoria/
    for folder in data/$categoria/*/; do
        # Pega o nome do ID (ex: SRR12345) removendo a barra final
        id=$(basename "$folder")
        
        echo "--- Processando QC para: $id ($categoria) ---"
        
        # Define os caminhos de entrada (procura os arquivos _1 e _2)
        R1="${folder}${id}_1.fastq"
        R2="${folder}${id}_2.fastq"
        
        # Define a pasta de saída específica para esta amostra
        OUT_SAMPLE="results/qc/$categoria/$id"
        mkdir -p "$OUT_SAMPLE/fastp"
        mkdir -p "$OUT_SAMPLE/fastqc"

        # 1. Rodar fastp
        fastp \
          -i "$R1" \
          -I "$R2" \
          -o "$OUT_SAMPLE/fastp/clean_${id}_R1.fastq" \
          -O "$OUT_SAMPLE/fastp/clean_${id}_R2.fastq" \
          -h "$OUT_SAMPLE/fastp/${id}_fastp.html" \
          -j "$OUT_SAMPLE/fastp/${id}_fastp.json"

        # 2. Rodar FastQC nos arquivos limpos
        fastqc \
          "$OUT_SAMPLE/fastp/clean_${id}_R1.fastq" \
          "$OUT_SAMPLE/fastp/clean_${id}_R2.fastq" \
          -o "$OUT_SAMPLE/fastqc"
    done
}

# Executar para as duas categorias
process_category "antigo"
process_category "recente"

# 3. Rodar MultiQC final unificando TUDO
echo "--- Gerando Relatório MultiQC Final ---"
multiqc results/qc -o results/qc/multiqc_final

echo "Pipeline de QC em lote finalizado!"

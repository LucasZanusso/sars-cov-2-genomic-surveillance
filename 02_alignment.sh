#!/bin/bash

# Criar pastas de saída
mkdir -p results/alignment/antigo
mkdir -p results/alignment/recente
REF="ref/wuhan_ref.fasta"

process_alignment() {
    local categoria=$1
    
    # Percorre as pastas onde o FASTP salvou os arquivos limpos
    for folder in results/qc/$categoria/*/; do
        id=$(basename "$folder")
        echo "--- Alinhando amostra: $id ($categoria) ---"

        R1="${folder}fastp/clean_${id}_R1.fastq"
        R2="${folder}fastp/clean_${id}_R2.fastq"
        OUT_DIR="results/alignment/$categoria/$id"
        mkdir -p "$OUT_DIR"

        # 1. Alinhamento com BWA MEM
        # -t 4 usa 4 threads (ajuste conforme seu PC)
        # -R adiciona o Read Group (essencial para alguns Variant Callers e boas práticas)
        bwa mem -t 4 -R "@RG\tID:$id\tSM:$id\tPL:ILLUMINA" "$REF" "$R1" "$R2" > "$OUT_DIR/${id}.sam"

        # 2. Converter SAM para BAM (binário)
        samtools view -S -b "$OUT_DIR/${id}.sam" > "$OUT_DIR/${id}.bam"

        # 3. Ordenar o BAM (Obrigatório para Variant Calling)
        samtools sort "$OUT_DIR/${id}.bam" -o "$OUT_DIR/${id}.sorted.bam"

        # 4. Indexar o BAM (Gera o arquivo .bai)
        samtools index "$OUT_DIR/${id}.sorted.bam"

        # 5. Limpeza: Remover arquivos intermediários pesados
        rm "$OUT_DIR/${id}.sam" "$OUT_DIR/${id}.bam"
        
        echo "Finalizado: $id"
    done
}

process_alignment "antigo"
process_alignment "recente"

echo "Todos os alinhamentos foram concluídos!"

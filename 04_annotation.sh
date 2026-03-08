#!/bin/bash

# Caminhos - AJUSTE O CAMINHO PARA O SEU snpEff.jar
SNPEFF_JAR="/home/lucas/Downloads/snpEff_5.4c/snpEff/snpEff.jar"
DB="NC_045512.2"

mkdir -p results/annotation/antigo
mkdir -p results/annotation/recente

process_annotation() {
    local categoria=$1
    
    # Percorre as pastas onde estão os VCFs filtrados
    for folder in results/vcf/$categoria/*/; do
        id=$(basename "$folder")
        echo "--- Anotando variantes: $id ($categoria) ---"

        VCF_IN="${folder}${id}_filtered.vcf"
        OUT_DIR="results/annotation/$categoria/$id"
        mkdir -p "$OUT_DIR"

        # Execução do snpEff
        # -v: verbose (mostra o que está fazendo)
        # -stats: gera o relatório HTML de estatísticas
        java -Xmx4g -jar "$SNPEFF_JAR" "$DB" "$VCF_IN" > "$OUT_DIR/${id}_annotated.vcf" \
             -stats "$OUT_DIR/${id}_summary.html"

        echo "Concluído: $id. Relatório gerado em ${id}_summary.html"
    done
}

process_annotation "antigo"
process_annotation "recente"

echo "Anotação funcional finalizada!"

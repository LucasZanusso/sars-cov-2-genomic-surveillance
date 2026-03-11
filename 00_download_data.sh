#!/bin/bash
set -euo pipefail

mkdir -p data/antigo data/recente logs

echo "--- Iniciando downloads ---"

while IFS=$'\t' read -r id grupo; do
    echo "Baixando $id em data/$grupo/$id..."
    mkdir -p "data/$grupo/$id"
    fasterq-dump "$id" --outdir "data/$grupo/$id" --split-files --progress >> logs/download.log 2>&1
done < samples.txt

echo "Download concluído! Estrutura de pastas organizada."

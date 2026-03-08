#!/bin/bash

# Criar estrutura básica
mkdir -p data/antigo data/recente logs

# Listas de IDs
ANTIGOS=("SRR24510773" "SRR24510774" "SRR24510775" "SRR24510776" "SRR24510777" "SRR24510778" "SRR24510780" "SRR24510779" "SRR24510781" "SRR24510782")
RECENTES=("SRR37461704" "SRR37461705" "SRR37461706" "SRR37461707" "SRR37461708" "SRR37461709" "SRR37461710" "SRR37461711" "SRR37461712" "SRR37461713")

echo "--- Iniciando download do Surto Antigo (2023) ---"
for id in "${ANTIGOS[@]}"; do
    echo "Baixando $id em data/antigo/$id..."
    # Cria a subpasta específica para o ID
    mkdir -p "data/antigo/$id"
    # Baixa os arquivos diretamente para dentro da subpasta
    fasterq-dump "$id" --outdir "data/antigo/$id" --split-files --progress >> logs/download.log 2>&1
done

echo "--- Iniciando download do Surto Recente (2026) ---"
for id in "${RECENTES[@]}"; do
    echo "Baixando $id em data/recente/$id..."
    # Cria a subpasta específica para o ID
    mkdir -p "data/recente/$id"
    # Baixa os arquivos diretamente para dentro da subpasta
    fasterq-dump "$id" --outdir "data/recente/$id" --split-files --progress >> logs/download.log 2>&1
done

echo "Download concluído! Estrutura de pastas organizada."

#!/bin/bash

# Criar pastas para salvar os arquivos VCF (Variant Call Format)
mkdir -p results/vcf/antigo
mkdir -p results/vcf/recente
REF="ref/wuhan_ref.fasta"

process_variants() {
    local categoria=$1
    
    # Percorre as pastas de alinhamento
    for folder in results/alignment/$categoria/*/; do
        id=$(basename "$folder")
        echo "--- Buscando mutações na amostra: $id ($categoria) ---"

        # Arquivo de entrada (BAM ordenado)
        BAM="${folder}${id}.sorted.bam"
        
        # Pasta de saída para o VCF
        OUT_DIR="results/vcf/$categoria/$id"
        mkdir -p "$OUT_DIR"

        # Passo 1: Empilhar as leituras (mpileup) e Chamar as Variantes (call)
        # -a DP: Garante que a Profundidade (Depth) seja anotada no arquivo para usarmos no filtro
        # -m: Permite múltiplas alternativas de alelos
        # -v: Mostra apenas as linhas onde há variação (ignora o que é igual à referência)
        bcftools mpileup -a DP -f "$REF" "$BAM" | bcftools call --ploidy 1 -mv -Ob -o "$OUT_DIR/${id}_raw.bcf"

        # Passo 2: Aplicar os Filtros do Usuário (QUAL > 30 e DP > 10)
        # -i 'QUAL>30 && DP>10': Inclui apenas variantes que passam nesta regra
        # -O v: Saída no formato VCF legível por humanos (texto padrão)
        bcftools filter -i 'QUAL>20 && INFO/DP>3' -O v -o "$OUT_DIR/${id}_filtered.vcf" "$OUT_DIR/${id}_raw.bcf"

        # (Opcional) Passo 3: Criar um arquivo indexado .gz para visualização em softwares como o IGV
        bgzip -c "$OUT_DIR/${id}_filtered.vcf" > "$OUT_DIR/${id}_filtered.vcf.gz"
        tabix -p vcf "$OUT_DIR/${id}_filtered.vcf.gz"

        # Limpar arquivo bruto para economizar espaço
        rm "$OUT_DIR/${id}_raw.bcf"

        echo "Variantes filtradas salvas em: ${id}_filtered.vcf"
    done
}

# Rodar para as duas categorias
process_variants "antigo"
process_variants "recente"

echo "Variant Calling finalizado com sucesso!"

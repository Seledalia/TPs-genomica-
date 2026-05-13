#!/usr/bin/env bash
# ============================================================
# Trabajo Práctico - Ensamblaje Genómico
# Adaptado para GitHub Codespace (Ubuntu)
# ============================================================
set -euo pipefail   # detiene el script si hay error

# ============================================================
# 0. PREPARACIÓN DEL AMBIENTE
# ============================================================
echo ">>> [0] Instalando micromamba..."

# Descargar e instalar micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
export PATH="$(pwd)/bin:$PATH"

# Crear entorno conda con todas las herramientas necesarias
./bin/micromamba create -y -n assembly \
    -c conda-forge -c bioconda \
    fastqc seqkit cutadapt minia spades taxonkit blast goldrush

# Instalar sra-tools desde apt
sudo apt-get update -q
sudo apt-get install -y -q sra-toolkit

# Configurar caché de SRA
mkdir -p ~/.ncbi
WORKDIR="$(pwd)"
micromamba run -n assembly vdb-config \
    -s /repository/user/main/public/root="${WORKDIR}/sra_cache"
mkdir -p sra_cache

echo ">>> Ambiente listo."

# ============================================================
# 1. OBTENCIÓN Y LIMPIEZA DE LECTURAS CORTAS (Illumina)
# ============================================================
echo ""
echo ">>> [1] Descargando lecturas de SRA (SRR10971381)..."

mkdir -p reads
NUM_READS=100000

fastq-dump -X "${NUM_READS}" \
    -A SRR10971381 \
    --split-files \
    --gzip \
    --outdir reads

echo "Archivos descargados:"
ls reads/

# Rutas a los archivos de lecturas
F1="reads/SRR10971381_1.fastq.gz"
F2="reads/SRR10971381_2.fastq.gz"

# Control de calidad inicial
echo ""
echo ">>> Stats de lecturas crudas:"
micromamba run -n assembly seqkit stats "${F1}" "${F2}"

# FastQC sobre lecturas crudas
echo ""
echo ">>> Corriendo FastQC sobre lecturas crudas..."
micromamba run -n assembly fastqc "${F1}" "${F2}" --outdir reads/

# ----- Trimming con Cutadapt -----
echo ""
echo ">>> Removiendo adaptadores y lecturas de baja calidad (cutadapt)..."

F1TRIM="reads/SRR10971381_trimmed_1.fastq.gz"
F2TRIM="reads/SRR10971381_trimmed_2.fastq.gz"

ADAPTER_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
MINLEN=50
QUAL=20

micromamba run -n assembly cutadapt \
    -a "${ADAPTER_FWD}" \
    -A "${ADAPTER_REV}" \
    -o "${F1TRIM}" \
    -p "${F2TRIM}" \
    "${F1}" "${F2}" \
    --minimum-length "${MINLEN}" \
    -q "${QUAL}"

echo ""
echo ">>> Stats de lecturas recortadas:"
micromamba run -n assembly seqkit stats "${F1TRIM}" "${F2TRIM}"

# ============================================================
# 2. ENSAMBLAJE (Minia + SPAdes)
# ============================================================
mkdir -p resultados/minia resultados/spades

# ----- Minia -----
echo ""
echo ">>> [2a] Ensamblando con Minia..."

ls reads/*trimmed* > read_list

MINIA_OUT="resultados/minia/minia"

micromamba run -n assembly minia \
    -in read_list \
    -out "${MINIA_OUT}" \
    -kmer-size 75

echo "Archivos generados por Minia:"
ls resultados/minia/

# ----- SPAdes -----
echo ""
echo ">>> [2b] Ensamblando con SPAdes (esto puede demorar ~15 min)..."

SPADES_DIR="resultados/spades"

micromamba run -n assembly spades.py \
    -1 "${F1TRIM}" \
    -2 "${F2TRIM}" \
    -o "${SPADES_DIR}" \
    -k 21,33,55,75

echo "Archivos generados por SPAdes:"
ls "${SPADES_DIR}/"

# ============================================================
# 3. EVALUACIÓN DE CALIDAD DE LOS ENSAMBLAJES
# ============================================================
echo ""
echo ">>> [3] Evaluando calidad de los ensamblajes..."

MINIA_RES="resultados/minia/minia.contigs.fa"
SPADES_RES="${SPADES_DIR}/scaffolds.fasta"

micromamba run -n assembly seqkit stats "${MINIA_RES}" "${SPADES_RES}" --all

# Liberar espacio eliminando las lecturas crudas
echo ""
echo ">>> Liberando espacio (eliminando reads/)..."
rm -rf reads/

# ============================================================
# 4. ENSAMBLAJE DE LECTURAS LARGAS (Nanopore con Goldrush)
# ============================================================
echo ""
echo ">>> [4] Descargando lecturas largas (SRR31637353 - Nanopore)..."

mkdir -p long_reads
NUM_READS_LONG=100000

fastq-dump -X "${NUM_READS_LONG}" \
    -A SRR31637353 \
    --gzip \
    --outdir long_reads

echo ""
echo ">>> Ensamblando con Goldrush..."

GS=100000   # tamaño genómico esperado (~100 Kb)
THREADS=4   # ajustar según los recursos del Codespace

cd long_reads/
gunzip -f SRR31637353.fastq.gz

micromamba run -n assembly goldrush path-polish \
    reads=SRR31637353 \
    G="${GS}" \
    t="${THREADS}"

echo "Archivos generados por Goldrush:"
ls -lh

cd ..

# ============================================================
# 5. IDENTIFICACIÓN TAXONÓMICA (BLAST contra virus de referencia)
# ============================================================
echo ""
echo ">>> [5] Descargando base de datos de virus para BLAST..."

micromamba run -n assembly update_blastdb.pl \
    --decompress ref_viruses_rep_genomes

# Taxonomía
echo ">>> Descargando datos de taxonomía..."
wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz

mkdir -p ~/.taxonkit
cp names.dmp nodes.dmp delnodes.dmp merged.dmp ~/.taxonkit/

# BLAST local
echo ">>> Corriendo BLAST local..."
micromamba run -n assembly blastn \
    -db ref_viruses_rep_genomes \
    -query "${SPADES_RES}" \
    -outfmt 6 > blast.txt

# Mejor hit por secuencia
sort -k1,1 -k12,12g blast.txt | awk '!seen[$1]++' > best_hits_per_contig.tsv

# Taxonomía de los hits
cut -f2 best_hits_per_contig.tsv > acc.list

micromamba run -n assembly blastdbcmd \
    -db ref_viruses_rep_genomes \
    -entry_batch acc.list \
    -outfmt "%a %T" | awk '{print $2}' > acc_taxid.tsv

micromamba run -n assembly taxonkit lineage acc_taxid.tsv \
    | micromamba run -n assembly taxonkit reformat \
        -r "Unclassified" \
        -f "{f}" > results.taxonomy

echo ""
echo ">>> Resultados de taxonomía guardados en: results.taxonomy"
echo ""
echo "============================================================"
echo " TP finalizado. Resultados en:"
echo "   - resultados/minia/   -> ensamblaje Minia"
echo "   - resultados/spades/  -> ensamblaje SPAdes"
echo "   - long_reads/         -> ensamblaje Goldrush"
echo "   - results.taxonomy    -> clasificación taxonómica"
echo "============================================================"

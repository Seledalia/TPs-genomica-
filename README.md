# TP Ensamblaje Genómico — GitHub Codespace

Trabajo Práctico de ensamblaje genómico adaptado para correr en **GitHub Codespaces** (Ubuntu).

---

## Contenido del repositorio

```
assembly.sh       <- script principal con todos los pasos del TP
README.md         <- este archivo
```

---

## Cómo crear y usar el Codespace

### 1. Crear el Codespace

1. Entrá a tu repositorio en GitHub.
2. Hacé clic en el botón verde **`<> Code`**.
3. Elegí la pestaña **Codespaces** → `Create codespace on main`.
4. Esperá a que cargue el entorno (puede tardar ~1 min la primera vez).

### 2. Darle permisos de ejecución al script

Una vez dentro del Codespace, en la terminal escribí:

```bash
chmod +x assembly.sh
```

### 3. Correr el script completo

```bash
./assembly.sh
```

> ⚠️ El script demora bastante porque descarga datos de SRA y corre SPAdes (~15-20 min en total). Podés dejarlo corriendo.

### 4. Ver los resultados

Los resultados quedan en estas carpetas:

| Carpeta / Archivo     | Contenido                              |
|-----------------------|----------------------------------------|
| `resultados/minia/`   | Ensamblaje con Minia                   |
| `resultados/spades/`  | Ensamblaje con SPAdes                  |
| `long_reads/`         | Ensamblaje con Goldrush (Nanopore)     |
| `results.taxonomy`    | Clasificación taxonómica (BLAST)       |
| `best_hits_per_contig.tsv` | Mejor hit BLAST por contig        |

---

## Cómo guardar los cambios en GitHub

Después de correr el TP o editar archivos, guardá todo en GitHub con estos comandos:

```bash
# 1. Agregar todos los archivos nuevos/modificados
git add .

# 2. Hacer el commit con un mensaje descriptivo
git commit -m "Agrego resultados del TP de ensamblaje"

# 3. Subir los cambios a GitHub
git push
```

> 💡 Si solo querés guardar el script (sin los resultados pesados), podés hacer:
> ```bash
> git add assembly.sh README.md
> git commit -m "Agrego script del TP"
> git push
> ```

---

## Notas importantes

- El script usa `micromamba` para manejar el ambiente con todas las herramientas bioinformáticas. Se instala automáticamente al correr `assembly.sh`.
- El Codespace por defecto tiene **4 cores y 8 GB de RAM** — suficiente para este TP.
- Si el Codespace se desconecta a mitad, podés volver a abrirlo desde GitHub y los archivos seguirán ahí (mientras no lo elimines).
- Los archivos `.fastq.gz` son pesados; si no necesitás subirlos a GitHub, podés agregarlos al `.gitignore`.

### .gitignore recomendado

Creá un archivo `.gitignore` con este contenido para no subir archivos de datos pesados:

```
reads/
sra_cache/
long_reads/*.fastq
*.tar.gz
*.fastq
*.fastq.gz
ref_viruses_rep_genomes*
taxdump*
*.dmp
```

---

## Preguntas del TP

### Parte 1 — Lecturas cortas

1. ¿Por qué las lecturas de Illumina son de baja calidad inicialmente y cómo ayuda el trimming?
2. ¿Qué diferencias ves entre el ensamblaje de Minia y el de SPAdes en términos de N50, número de contigs y longitud total?

### Parte 2 — Lecturas largas

1. ¿Por qué con el mismo número de lecturas Nanopore se obtiene un ensamblaje de mejor calidad que con Illumina?
2. ¿Por qué Goldrush no requiere definir un tamaño de k-mer?
3. ¿Qué pasos adicionales ofrece Goldrush según su [repositorio en GitHub](https://github.com/BirolLab/goldrush)?

### Parte 3 — Identificación taxonómica

1. ¿Cuáles fueron los principales microorganismos en la muestra?
2. Modificá el script de Python para contar hits a nivel de especie u orden. ¿Hay algún grupo relevante?
3. ¿Los resultados son consistentes con una infección por coronavirus?
4. ¿Qué sesgos tiene este análisis? ¿Podrías diagnosticar solo con estos resultados?

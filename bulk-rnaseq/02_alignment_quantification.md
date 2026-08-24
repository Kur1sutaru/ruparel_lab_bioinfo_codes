# 02 — Alignment / Quantification and Building a Count Matrix

Two common paths are shown — pick one. **Salmon** (pseudo-alignment) is faster and the more common modern default for a straightforward DEG workflow. **STAR + featureCounts** gives you a real BAM if you need it for other purposes.

## Option A — Salmon (recommended default)

### Step 1 — Build a Salmon index (once per reference)

```bash
# Download a transcriptome FASTA + GTF for your species (e.g. from GENCODE/Ensembl) first
salmon index \
  -t transcriptome.fa.gz \
  -i salmon_index \
  -k 31 \
  -p 8
```

### Step 2 — Quantify each sample

```bash
mkdir -p salmon_quant
for r1 in trimmed/*_R1.trimmed.fastq.gz; do
  r2=${r1/_R1/_R2}
  sample=$(basename "$r1" _R1.trimmed.fastq.gz)
  salmon quant \
    -i salmon_index \
    -l A \
    -1 "$r1" -2 "$r2" \
    -p 8 \
    --validateMappings \
    --gcBias \
    -o "salmon_quant/${sample}"
done
```

`-l A` auto-detects library strandedness. `--gcBias` corrects for GC-content-driven quantification bias — generally worth including.

### Step 3 — Import into R with tximport

```r
library(tximport)
library(tidyverse)

samples <- list.dirs("salmon_quant", full.names = FALSE, recursive = FALSE)
files <- file.path("salmon_quant", samples, "quant.sf")
names(files) <- samples

# tx2gene: a 2-column data.frame/table mapping transcript ID -> gene ID,
# built from your GTF (e.g. with the `GenomicFeatures` package, or a pre-made
# tx2gene table from your reference provider)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

counts <- round(txi$counts)  # gene-level count matrix, samples as columns
```

## Option B — STAR + featureCounts

### Step 1 — Build a STAR genome index (once per reference)

```bash
STAR --runMode genomeGenerate \
  --genomeDir star_index \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile annotation.gtf \
  --sjdbOverhang 100 \
  --runThreadN 8
```

### Step 2 — Align each sample

```bash
mkdir -p star_out
for r1 in trimmed/*_R1.trimmed.fastq.gz; do
  r2=${r1/_R1/_R2}
  sample=$(basename "$r1" _R1.trimmed.fastq.gz)
  STAR --runThreadN 8 \
    --genomeDir star_index \
    --readFilesIn "$r1" "$r2" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "star_out/${sample}_"
  samtools index "star_out/${sample}_Aligned.sortedByCoord.out.bam"
done
```

### Step 3 — Count reads per gene with featureCounts

```bash
featureCounts \
  -a annotation.gtf \
  -o counts.txt \
  -T 8 -p --countReadPairs \
  star_out/*_Aligned.sortedByCoord.out.bam
```

### Step 4 — Load into R

```r
counts_raw <- read.delim("counts.txt", comment.char = "#")
counts <- counts_raw[, -(1:6)]        # drop Chr/Start/End/Strand/Length columns
rownames(counts) <- counts_raw$Geneid
colnames(counts) <- gsub("star_out/|_Aligned.*", "", colnames(counts))
```

## Checkpoint

A gene (rows) × sample (columns) integer count matrix in R, with sample names as column names matching your sample metadata table exactly — you'll need this alignment for Tutorial 03.

## Common pitfalls

- Mismatched sample names between the count matrix columns and your metadata table — DESeq2 will error or, worse, silently mismatch samples to conditions if column order is assumed rather than checked.
- Using an unstranded library type setting (`-l`/`--fr-unstranded` etc.) when your library prep was actually stranded, or vice versa — check your library prep kit documentation, or infer it from a tool like `RSeQC`'s `infer_experiment.py` if unsure.
- Forgetting `--countReadPairs`/paired mode flags for paired-end data in featureCounts, which silently produces wrong counts.

Next: [`03_deseq2_degs.md`](03_deseq2_degs.md)

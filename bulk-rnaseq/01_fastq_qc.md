# 01 — FASTQ Quality Control and Trimming

## Goal
Assess raw read quality and produce trimmed, analysis-ready FASTQs.

## Step 1 — Run FastQC on all raw FASTQs

```bash
mkdir -p qc/fastqc_raw
fastqc raw_fastq/*.fastq.gz -o qc/fastqc_raw -t 4
```

## Step 2 — Aggregate with MultiQC

```bash
multiqc qc/fastqc_raw -o qc/fastqc_raw
```

Open `qc/fastqc_raw/multiqc_report.html`. What to actually look at:
- **Per-base sequence quality** — should stay comfortably in the green (Phred ≥ 28-30ish) across most of the read length; a quality drop at the read's 3' end is normal and handled by trimming.
- **Adapter content** — nonzero adapter contamination is expected and fine; trimming (Step 3) removes it.
- **Per-sequence GC content** — should roughly match a single peak near your organism's expected GC%; a bimodal or shifted distribution can indicate contamination (e.g. rRNA, another organism).
- **Overrepresented sequences** — check whether flagged sequences are just adapters (fine) or something unexpected (rRNA, a specific transcript suggesting degradation).

## Step 3 — Trim adapters and low-quality bases with fastp

```bash
mkdir -p trimmed qc/fastp

for r1 in raw_fastq/*_R1.fastq.gz; do
  r2=${r1/_R1/_R2}
  sample=$(basename "$r1" _R1.fastq.gz)
  fastp \
    -i "$r1" -I "$r2" \
    -o "trimmed/${sample}_R1.trimmed.fastq.gz" \
    -O "trimmed/${sample}_R2.trimmed.fastq.gz" \
    --detect_adapter_for_pe \
    --json "qc/fastp/${sample}.json" \
    --html "qc/fastp/${sample}.html" \
    --thread 4
done
```

`fastp` autodetects adapters for paired-end data and handles quality trimming in one pass; adjust `-i`/`-I` (single- vs paired-end) to match your library type.

## Step 4 — Re-run QC on trimmed reads

```bash
mkdir -p qc/fastqc_trimmed
fastqc trimmed/*.fastq.gz -o qc/fastqc_trimmed -t 4
multiqc qc/fastqc_trimmed -o qc/fastqc_trimmed
```

Confirm adapter content has dropped to near-zero and quality scores look clean across the full read length.

## Checkpoint

A `trimmed/` folder with paired trimmed FASTQs per sample, and a MultiQC report showing clean per-base quality with adapter content removed. Note the % reads surviving trimming per sample (fastp reports this) — a sample losing an unusually large fraction of reads is worth flagging before proceeding.

## Common pitfalls

- Trimming too aggressively (e.g. very strict quality cutoffs) and losing usable read length — fastp's defaults are reasonable for most standard Illumina data; don't over-tune without a reason.
- Mismatching R1/R2 file pairs when looping over samples — always verify pairing (e.g. by read count) before proceeding to alignment.
- Skipping the post-trim QC step and assuming trimming worked — always re-check.

Next: [`02_alignment_quantification.md`](02_alignment_quantification.md)

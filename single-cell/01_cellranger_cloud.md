# 01 — Running Cell Ranger on 10x Genomics Cloud Analysis

## Goal
Go from raw FASTQ files to a feature-barcode matrix, entirely in the browser, with no local Cell Ranger install.

## Concepts

10x Genomics Cloud Analysis runs the standard `cellranger count` pipeline on 10x's infrastructure. It demultiplexes reads, aligns to a reference transcriptome (STAR under the hood), calls real cells vs. empty droplets, and outputs a feature-barcode matrix plus a QC summary (`web_summary.html`) — the same outputs you'd get from a local Cell Ranger install, without needing a workstation with 64+ GB RAM.

## Step 1 — Create an account

1. Go to [https://www.10xgenomics.com/analysis-cloud](https://www.10xgenomics.com/analysis-cloud) (or [cloud.10xgenomics.com](https://cloud.10xgenomics.com) if you already know the direct link).
2. Click **Sign Up** and register with your institutional email.
3. Verify your email and log in. You'll land on the Cloud Analysis dashboard.
4. Every institution/lab typically shares one **organization** — ask whoever set up the lab's account to add you as a member so uploaded data and run history are shared, rather than creating a personal silo.

## Step 2 — Know your inputs before uploading

You need, per sample:
- **FASTQ files** following 10x naming convention: `SampleName_S1_L00X_R1_001.fastq.gz`, `..._R2_001.fastq.gz` (and `_I1_`/`_I2_` if present). If your sequencing core gave you differently-named files, use `bcl2fastq`/`cellranger mkfastq` conventions or ask the core to re-demux — Cloud Analysis expects this naming to auto-detect samples.
- **Reference genome/transcriptome** — 10x provides prebuilt references (human GRCh38, mouse mm10, etc.) directly in the Cloud interface picker. For non-standard species (e.g. marmoset, rat with a custom annotation), you'll need to build a custom reference locally with `cellranger mkref` and upload it, or check if 10x has since added a prebuilt one for your species.
- **Chemistry version** (v2 / v3 / Flex / etc.) — check your library prep kit; Cloud Analysis usually auto-detects this but verify it against what the core reports.

## Step 3 — Upload FASTQs

1. From the dashboard, click **Upload Data** (or **New Run** → it will prompt you to upload first if you haven't already).
2. Drag and drop the FASTQ files, or connect a cloud bucket (S3/GCS) if your core deposits sequencing data there directly — this avoids downloading large files locally just to re-upload them.
3. Wait for upload + checksum verification to complete. Large runs (multiple lanes, deep sequencing) can take a while — this is a good stopping point.

## Step 4 — Configure and launch the run

1. Click **New Analysis** → choose **Cell Ranger** → **Gene Expression (count)** (or **Multi** if you have multiplexed/hashed/CITE-seq samples — see note below).
2. Select your uploaded FASTQs for this sample.
3. Select the reference (prebuilt or your custom uploaded one).
4. Leave **expected cell count** at the pipeline default unless you have a strong prior from your loading concentration — over-specifying this can bias cell calling.
5. Name the run something you'll recognize later (e.g. `tongue_human_rep1_20260824`) — you will thank yourself during Tutorial 02 when loading a dozen samples.
6. Click **Launch**.

> **Multiplexed / hashed samples:** if multiple biological samples were pooled into one 10x lane (cell hashing, CMO multiplexing, or Flex), use the **Cell Ranger Multi** pipeline instead of plain `count`, and provide the hashtag/CMO-to-sample mapping in the config step. This demultiplexes samples back apart as part of the same run.

## Step 5 — Monitor and retrieve outputs

- Runs typically take a few hours depending on read depth and cell count; you'll get an email/dashboard notification on completion.
- Open the run and download (or browse in-browser):
  - `web_summary.html` — **read this first**. Check: estimated number of cells, median genes/UMIs per cell, sequencing saturation, fraction reads in cells. A low fraction-reads-in-cells or unusual cell number vs. expected loading is worth flagging before moving on.
  - `filtered_feature_bc_matrix/` (or `.h5`) — the matrix you'll load into Seurat in Tutorial 02. Use the **filtered** (not raw) matrix unless you specifically need to inspect empty droplets yourself.
  - `metrics_summary.csv` — machine-readable version of the web summary, useful for building a QC table across many samples.

## Checkpoint

You should have, per sample, either a downloaded `filtered_feature_bc_matrix.h5` or a `filtered_feature_bc_matrix/` folder (containing `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`), plus a `web_summary.html` you've actually opened and looked at.

## Common pitfalls

- **FASTQ naming mismatch** — Cloud Analysis silently fails to group files into one sample if names deviate from the 10x convention. Rename before uploading.
- **Wrong chemistry auto-detected** — always double check the detected chemistry in the run config against what your core reports.
- **Using the raw instead of filtered matrix** downstream — the raw matrix includes hundreds of thousands of empty-droplet barcodes and is not what you load into Seurat for a standard analysis.
- **Not reading `web_summary.html`** before moving on — catching a bad run (e.g. very low cells detected, high ambient RNA) here saves hours of debugging in Tutorial 02.

Next: [`02_qc_filtering.md`](02_qc_filtering.md)

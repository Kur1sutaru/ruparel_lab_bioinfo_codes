We will use our lab published bulk-RNA-seq data for this practice, available at 

Ibrahim, T., Wu, P., Wang, LJ. et al. Sex-dependent differences in the genomic profile of lingual sensory neurons in naïve and tongue-tumor bearing mice. Sci Rep 13, 13117 (2023). https://doi.org/10.1038/s41598-023-40380-6


These are RNA-seq of FACS-sorted, WGA-retro-labeled tongue-innervating trigeminal ganglion neurons from Balb/c nude mice, 4 groups (Male Normal, Male Tumor, Female Normal, Female Tumor), n=3 each, SMART-seq2 libraries, 50bp single-end reads. I've built the pipeline around that design — 12 samples, sex × tumor-status comparisons (FN vs MN, MT vs MN, FT vs FN, and FT vs MT if you want it).

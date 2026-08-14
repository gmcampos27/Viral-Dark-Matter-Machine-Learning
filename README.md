# Viral-Dark-Matter-Machine-Learning
Viral Dark Matter classification model

This project implements a pipeline for the classification of viral "Dark Matter" in metagenomic datasets.
It combines Genomic Features with Machine Learning, we classify sequences into 34 viral families, providing a fast and reliable alternative to traditional alignment-based methods.

In metagenomics, a significant portion of sequences remains unclassified (Dark Matter).
This repository was specifically designed to explore viral diversity in plasma samples, focusing on Anelloviridae.

Features: Extraction of GC content, Tetranucleotide Frequencies (4-mers), and Pfam functional domains.
Dual Modeling: Comparative analysis between Random Forest and XGBoost.
Statistics Analyses: Model validation using 1000x Bootstrap, ANOVA, and Cohen’s $d$ effect size.
Taxonomic Validation: Integrated R scripts for ggtree phylogenetics and ICTV-compliant pairwise identity analysis (69% threshold).

├── pipeline/            # Quality control, host depletion & assembly
├── codes_ml/            # Feature extraction, training (RF/XGB) & stats
├── Anellovirus_Tree/    # Phylogeny (ggtree) & Pairwise Identity plot
├── Models/              # Trained .pkl models & LabelEncoders

Please go to https://doi.org/10.5281/zenodo.21907977 to see the dataset

Contigs are avaliable at the following Accession Numbers: GenBank PZ533894-PZ533937 

Public Avaliable Samples utilized in this study: 
SRR32222465
SRR32222464
SRR32222453
SRR32222442
SRR32222431
SRR32222420
SRR32222409
SRR32222408
SRR32222407
SRR32222406
SRR32222463
SRR32222462
SRR32222461
SRR32222460
SRR32222459

Under the Bioproject PRJNA1219095 https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1219095/

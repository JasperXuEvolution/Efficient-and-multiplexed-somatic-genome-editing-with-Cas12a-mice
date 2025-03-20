# Efficient and multiplexed somatic genome editing with Cas12a mice
This repository contains data analysis scripts for the study "Efficient and Multiplexed Somatic Genome Editing with Cas12a Mice."

## Overview
Somatic genome editing in mouse models has increased our understanding of the in vivo effects of genetic alterations in areas ranging from neuroscience to cancer biology and beyond. However, existing models are limited in their ability to create multiple targeted edits. Thus, our understanding of the complex genetic interactions that underlie development, homeostasis, and disease remains incomplete. Cas12a is an RNA-guided endonuclease with unique attributes that enable simple targeting of multiple genes with crRNA arrays containing tandem guides. To accelerate and expand the generation of complex genotypes in somatic cells, we generated transgenic mice with Cre-regulated and constitutive expression of enhanced Acidaminococcus sp. Cas12a (enAsCas12a). In these mice, enAsCas12a-mediated somatic genome editing robustly generated compound genotypes, as exemplified by the initiation of diverse cancer types driven by homozygous inactivation of trios of tumor suppressor genes or an oncogenic translocation. We further integrated these modular crRNA arrays with clonal barcoding to quantify the size and number of tumors with each array, as well as the efficiency of each crRNA. Efficiency and TSG combos These Cas12a alleles will enable the rapid generation of disease models and broadly facilitate the high-throughput investigation of coincident genomic alterations in somatic cells in vivo.

## Contact
For questions or comments, please reach out to xhq@stanford.edu.

## Cancer model
- Scripts for analyzing Cas12a Ultra-seq data in three different tumor models:
  - Oncogene-negative tumors
  - Pancreatic ductal adenocarcinoma (PDAC)
  - Small cell lung cancer (SCLC)

- PE150 folder: Scripts for PE150 sequencing data (using SCLC as an example).
- PE300 folder: Scripts for PE300 sequencing data.

## Cas12a_Efficiency 
- Scripts for analyzing Cas12a Ultra-seq data from efficiency experiments.

## Triple_Knockout 
- Scripts for analyzing Cas12a Ultra-seq data from the triple tumor suppressor knockout experiment.

## Figure_Plot
- Scripts for generating all Cas12a Ultra-seq experiment-related figures.


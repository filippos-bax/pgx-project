A small end-to-end NGS to PGx analysis pipeline using real sequencing data

This project implements a complete, reproducible workflow for calling and annotating CYP2D6 variants from real Illumina amplicon sequencing data. It uses standard clinical and research bioinformatics tools (FastQC, BWA, SAMtools, GATK, bcftools, VEP) and performs gene-focused pharmacogenomic interpretation.

The purpose of the project is to demonstrate practical competence with NGS pipelines in a PGx context and to produce an annotated CYP2D6 variant report.

SRA data acquisition: SRA Toolkit - prefetch, fasterq-dump
Quality control: FastQC
Reference genome acquisition: chromosome selection w/ samtools faidx, indexing w/ bwa idx
Alignment: bwa-mem
Variant calling: GATK HaplotypeCaller
Variant filtering GATK VariantFiltration
Region extraction: bcftools
Functional annotation: VEP
VEP file parsing & visualization: Python

This pipeline identified the existence of clinically significant variants in the sequencing data used, namely rs3892097, rs1135840 and rs1065852, which are found in alleles of pharmakogenomic interest that are related to decreased or absent drug metabolism.
# Genotype Processing Scripts  

This repository provides **Bash** and **R** scripts to efficiently process and analyze genotype files.  

## Features  

### 1. Genotype File Processing (SNP_in_cols.sh and edit_genotypes.sh) 
- Work with genotype files, including `chromosome.data` and `imputed_genotypefiles`.   

### 2. File Comparison (edit_genotypes.sh)  
- Compare genotype files to:  
  - Identify **SNP removals**.  
  - Detect updates in genotype files.  

### 3. Quality Control (QC-SNP.R)  
Perform quality control checks on genotype data, including:  
- **Call Rate**: Ensure sufficient SNP call coverage.  
- **Call Frequency**: Verify the frequency of SNP calls.  
- **Minor Allele Frequency (MAF)**: Exclude SNPs with low-frequency alleles.  

### 4. Kernel Computation  
- Compute:  
  - **Linear kernels**.  
  - **Gaussian kernels**.  

## Requirements  

### Software  
- **Bash**: For processing pipeline scripts.  
- **R**: Version 4.0+ recommended.  

#!/bin/bash

# PART 1: Create list of SNP positions to be removed from genotype files

# Define file paths
old_SNPfile="/Data1/chromosome.data"
new_SNPfile="/Data2/chromosome.data"
output_remove_old="remove_from_old"
output_remove_new="remove_from_new"

# Step 1: Extract the SNP names (first column) from both files
awk '{print $1}' "$old_SNPfile" > old_snps
awk '{print $1}' "$new_SNPfile" > new_snps

# Step 2: Find common SNPs
sort old_snps new_snps | uniq -d > common_snps

# Step 3: Find SNPs to remove from each file
# SNPs in old file but not in common
grep -Fvxf common_snps old_snps > remove_new
# SNPs in new file but not in common
grep -Fvxf common_snps new_snps > remove_old

# Step 4: Generate output files sorted by Overall position in descending order
awk 'NR==FNR {a[$1]; next} $1 in a' remove_new "$old_SNPfile" | sort -k3,3nr > "$output_remove_old"
awk 'NR==FNR {a[$1]; next} $1 in a' remove_old "$new_SNPfile" | sort -k3,3nr > "$output_remove_new"

# Clean up intermediate files
rm old_snps new_snps common_snps remove_new remove_old

echo "Files created:"
echo "- $output_remove_old"
echo "- $output_remove_new"


# PART 2: Remove SNPs from genotype files using the list of removal positions

# SNP file 1 - Old SNP file updated due to potential marker remotions 
input="/Data1/imputed_genotypes1.txt"

# Count markers in original file
awk '{print $4}' "$input"  | awk '{print length($0)}' | head -n1

# Open the list of SNP to be remove from the SNP file
awk '{print $3}' remove_from_old > colsfile
awk '{print $4}' "$input" > snpfile

# Loop in SNP file to remove the cols in list
for col in $(cat colsfile)
 do

  cut -b $col --complement snpfile > tmpr
  mv tmpr snpfile

done

# Paste ids to the SNP file
awk '{print length($1)}' snpfile | head -n1
awk '{print $1}' "$input" > ids
paste -d ' ' ids snpfile > old_genotypes.txt

# Count markers in corrected file
awk '{print $2}' old_genotypes.txt | awk '{print length($0)}' | head -n1

# Clean up intermediate files
rm snpfile colsfile tmpr ids


# SNP file 2 - New SNP file adjusted to match with the old 
input="/Data2/imputed_genotypes2.txt"

# Count markers in original file
awk '{print $4}' "$input"  | awk '{print length($0)}' | head -n1

# Open the list of SNP to be remove from the SNP file
awk '{print $3}' remove_from_new > colsfile
awk '{print $4}' "$input" > snpfile

# Loop in SNP file to remove the cols in list
for col in $(cat colsfile)
 do

  cut -b $col --complement snpfile > tmpr
  mv tmpr snpfile

done

# Paste ids to the SNP file
awk '{print length($1)}' snpfile | head -n1
awk '{print $1}' "$input" > ids
paste -d ' ' ids snpfile > new_genotypes.txt

# Count markers in corrected file
awk '{print $2}' new_genotypes.txt | awk '{print length($0)}' | head -n1

# Clean up intermediate files
rm snpfile colsfile tmpr ids

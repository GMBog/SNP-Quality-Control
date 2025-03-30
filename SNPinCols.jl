# Count lines in the files
countlines("GenotypeFile.snp")

# Read the SNP file and process columns
genotype_data = readlines("GenotypeFile.snp");

# Extract first column (IDs)
ids = [split(line)[1] for line in genotype_data];

# Extract second column (SNPs) and insert spaces between characters
snps = [join(split(split(line)[2], ""), " ") for line in genotype_data];

# Create new SNP file
open("Genotypes.snp", "w") do io
    for (id, snp) in zip(ids, snps)
        println(io, id * " " * snp)
    end
end

# Get the max number of columns (=markers)
maximum(length.(split.(readlines("Genotypes.snp"))))

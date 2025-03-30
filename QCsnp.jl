import Pkg;
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Statistics")
Pkg.add("DelimitedFiles")

using CSV
using DataFrames
using Statistics
using DelimitedFiles

# Read SNP markers data
SNPdata = CSV.read("snpfile.dat", DataFrame; header = false);
size(SNPdata)

# Replace values 3, 4, and 5 with `missing`
for col in names(SNPdata)
    if eltype(SNPdata[!, col]) <: Number  # Apply only to numeric columns
        # Ensure column supports missing values
        SNPdata[!, col] = convert(Vector{Union{eltype(SNPdata[!, col]), Missing}}, SNPdata[!, col])

        # Replace values 3, 4, and 5 with `missing`
        SNPdata[!, col][in([3, 4, 5]).(SNPdata[!, col])] .= missing
    end
end

# SNP Call Rate (CR95%) 
# Check for missing frequency in each column and remove columns (SNP) based on the threshold 
MissSNP = Float64[]  # Initialize an empty array to store the results
for col in names(SNPdata)
    missing_proportion = sum(ismissing.(SNPdata[!, col])) / nrow(SNPdata)
    push!(MissSNP, missing_proportion)
end

sumCR = DataFrame(
           Column = names(SNPdata),
           Percentage_Missing = MissSNP .* 100
       )

SNPremove = MissSNP .>= 0.95;


# SNP Call Freq (CF95%)
# Check for missing frequency in each row and remove rows(Animal) based on the threshold 
MissANI = Float64[]  # Initialize an empty array to store the results
for row in 1:nrow(SNPdata)
    missing_proportion = sum(ismissing.(Vector(SNPdata[row, :]))) / ncol(SNPdata)
    push!(MissANI, missing_proportion)
end

sumCF = DataFrame(
           Row = SNPdata[:,1],
           Percentage_Missing = MissANI .* 100
       )

ANIremove = MissANI .>= 0.95;


# Removing SNPs and Animals from data
SNPdata = SNPdata[:,.!SNPremove]
SNPdata = SNPdata[.!ANIremove,:]
size(SNPdata)


# Replacing missings with the mean genotype per SNP 
SNP_impdata = copy(SNPdata)
for col in names(SNP_impdata)[2:end]  # Skip the first column (index 1)
    mean_value = mean(skipmissing(SNP_impdata[!, col]))
    SNP_impdata[ismissing.(SNP_impdata[!, col]), col] .= round(mean_value)
end


# Minor Allele Frequency (MAF)
# Checking for minor allele frequency (1%<x<99%)
MAF = Float64[]

for col in 2:ncol(SNP_impdata)
    freq = sum(SNP_impdata[:, col]) / (2 * nrow(SNP_impdata))
    push!(MAF, freq)  # Store the frequency in the MAF vector
end

snpremove = (MAF .<= 0.01) .| (MAF .>= 0.99);
snpremove = vcat(false, snpremove)

SNP_impdata = SNP_impdata[:, .!snpremove];
size(SNP_impdata)

# Order the SNP data based on the GenotypeID (first column)
SNP_impdata = sort!(SNP_impdata, by = x -> x[1]);
SNP_impdata[1:10,1:10]

# Save imputed SNP data to a file
writedlm("imputed_file.snp", Matrix(SNP_impdata), ' ')

# Computing genomic relationship matrix (G)

# Read genotypic data
X = Matrix(SNP_impdata[:, 2:end])

p = size(X, 2)

# Scale matrix to center in 0 and variance in 1
S = (X .- mean(X, dims=1)) ./ std(X, dims=1)

# Computing G matrix as the crossproduct between the scale matrix divided number of markers (S%*%t(S))/p
G = (S * S') / p
size(G)

# Save G matrix
writedlm("G_LinearK", G, ' ')

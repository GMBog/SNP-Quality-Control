using DataFrames
using CSV

function process_SNP_files(chrmap1_path::String, chrmap2_path::String, imputed_genotypes1_path::String, imputed_genotypes2_path::String)

    # Step 1: Read chromosome data files 
    chrmap1 = CSV.read(chrmap1_path, DataFrame)
    chrmap2 = CSV.read(chrmap2_path, DataFrame)

    # Step 2: Find common SNP names in chromosome data
    commonSNP = intersect(chrmap1.Name, chrmap2.Name)

    # Step 3: Create df with markers to remove from each SNP file
    remove_from_chrmap1 = chrmap1[in.(chrmap1.Name, Ref(setdiff(chrmap1.Name, commonSNP))), :]
    remove_from_chrmap1 = sort(remove_from_chrmap1, "Overall", rev = false)

    remove_from_chrmap2 = chrmap2[in.(chrmap2.Name, Ref(setdiff(chrmap2.Name, commonSNP))), :]
    remove_from_chrmap2 = sort(remove_from_chrmap2, "Overall", rev = false)

    # Process imputed genotypes file 1
    lines1 = readlines(imputed_genotypes1_path)

    function parse_line(line)
        m = match(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(.+)$", line)
        return isnothing(m) ? nothing : (parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]), strip(m[4]))
    end

    df1 = filter(!isnothing, parse_line.(lines1))
    df1 = DataFrame(:Key => first.(df1), :Col2 => getindex.(df1, 2), :Col3 => getindex.(df1, 3), :Col4 => last.(df1))

    snp_matrix1 = [collect(snp) for snp in df1[:, 4]]
    SNPinCols1 = permutedims(DataFrame(snp_matrix1, :auto))

    keep_indices1 = setdiff(1:size(SNPinCols1, 2), collect(remove_from_chrmap1.Overall))
    SNPinCols1 = SNPinCols1[:, keep_indices1]

    df1 = hcat(df1[:, 1:3], SNPinCols1)

    # Process imputed genotypes file 2
    lines2 = readlines(imputed_genotypes2_path)

    df2 = filter(!isnothing, parse_line.(lines2))
    df2 = DataFrame(:Key => first.(df2), :Col2 => getindex.(df2, 2), :Col3 => getindex.(df2, 3), :Col4 => last.(df2))

    snp_matrix2 = [collect(snp) for snp in df2[:, 4]]
    SNPinCols2 = permutedims(DataFrame(snp_matrix2, :auto))

    keep_indices2 = setdiff(1:size(SNPinCols2, 2), collect(remove_from_chrmap2.Overall))
    SNPinCols2 = SNPinCols2[:, keep_indices2]

    df2 = hcat(df2[:, 1:3], SNPinCols2)

    return df1, df2
end


#df1, df2 = process_SNP_files("data/chromosome1.data", "data/chromosome2.data", "data/imputed_genotypes1.txt", "data/imputed_genotypes2.txt");
#size(df1)
#size(df2)

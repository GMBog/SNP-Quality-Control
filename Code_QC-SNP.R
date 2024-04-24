### Quality control and computation of genomic and microbiome matrices
### Computation of Linear and Gaussian Kernels for BGLR

library(tidyr)
library(dplyr)

## SNP file should have in the first column animal ID and then SNP markers (one per column)
## To do that the SNP_in_columns.sh split SNP markers in columns!

## Read SNP markers 

SNPdata = read.table("~/SNPdata.txt", quote="\"", comment.char="")

### Imputation of missing values (3, 4 and 5) with NAs first and then with the mean of each column
SNPdata[SNPdata == 3] <- NA
SNPdata[SNPdata == 4] <- NA
SNPdata[SNPdata == 5] <- NA

### Check for NA frequency in each column and remove SNP with higher than 0.05 (SNP call rate 95%) 
MissSNP = colSums(is.na(SNPdata))/nrow(SNPdata)
table(MissSNP>=0.95)

SNPremove = MissSNP>=0.95
str(SNPremove)

### Check for NA frequency in each row and remove ANIMAL with higher than 0.05 (ANIMAL call rate 95%) 
MissANI = rowSums(is.na(SNPdata))/ncol(SNPdata)
table(MissANI>=0.95)

ANIremove = MissANI>=0.95
str(ANIremove)

### Removing SNP and ANIMALS
SNPdata = SNPdata[,!SNPremove]
SNPdata = SNPdata[!ANIremove,]
dim(SNPdata)

### Replace NAs with mean genotype per SNP 
SNP_impdata = SNPdata %>% mutate_all(~ifelse(is.na(.x), round(mean(.x, na.rm = TRUE),0), .x))

### Check for minor allele frequency
SNPMAF = SNP_impdata[,-1]
freq = colSums(SNPMAF)/(2*nrow(SNPMAF))
snpremove = freq <= 0.01 | freq >= 0.99
table(snpremove)

SNPMAF = SNPMAF[,!snpremove]
dim(SNPMAF)
SNP_impdata = cbind(SNP_impdata[1],SNPMAF)

### Order SNP data based on the GenotypeID
SNP_impdata = SNP_impdata[order(SNP_impdata[,1]),]
head(SNP_impdata[1])

### Save imputed SNP data
write.table(SNP_impdata, file="~/SNP_imputedQC.txt", row.names = F, col.names = F, quote = F)



### Computing G relationship matrix assuming LINEAR kernel

#Read genotypic data
X = as.matrix(SNP_impdata[,-1])

#Number of markers
p = ncol(X)
p

#Scale matrix to center in 0 and variance in 1
S = scale(X, center =  TRUE, scale = TRUE)

#Compute the G matrix as the crossproduct between the scale matrix divided number of markers (S%*%t(S))/p
G = tcrossprod(S)/p
dim(G)

#Save G matrix
write.table(G, file="~/G_LinearK", row.names = F, col.names = F, quote = F)


### Computing G relationship matrix assuming GAUSSIAN kernel

# Read genotypic data
X = as.matrix(SNP_impdata[,-1])

#Scale matrix to center in 0 and variance in 1
S = scale(X, center =  TRUE, scale = TRUE)

#Computing Euclidean distance between SNPs
D = as.matrix(dist(S,method="euclidean"))^2
D = D/mean(D)

#Save D matrix
write.table(D, file="~/DG", row.names = F, col.names = F, quote = F)

#Computing G matrix using a h-value of 0.7 
h = 1.5
GK = exp(-h*D)

#Save G matrix
write.table(GK, file="~/G_GaussianK", row.names = F, col.names = F, quote = F)



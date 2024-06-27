#!/bin/bash
#SBATCH --job-name=SNPdata
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50gb                    #memory per each job
#SBATCH --time=90:00:00                                 #time per each job

#Open the SNP file
awk '{print $1}' ~/dirfiles/GenotypesFile > Ids
awk '{print $2}' ~/dirfiles/GenotypeFile | sed 's/./& /g' > SNP_in_columns.txt

#Paste ids to the SNP file
paste -d ' ' Ids SNP_in_columns.txt > SNPdata.txt

rm SNP_in_columns.txt
rm Ids

module load R
Rscript QC-SNP.R

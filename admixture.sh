#!/bin/bash

# Convert the VCF file to BED format
plink --allow-extra-chr --vcf /home/roo/data/anlaysis_revise/dataset2/dataset2.recode.vcf --make-bed --out dataset2

#!/bin/bash

# Set the input file and output prefix
input_file="dataset2.bed"
output_prefix="output"
num_runs=100

# Loop over different K values from 2 to 10
for K in {2..10}; do
    # Loop over the number of runs
    for run in $(seq 1 $num_runs); do
        # Run admixture command and redirect output to log file
        admixture -j4 --cv "$input_file" "$K" > "${output_prefix}_K${K}_run${run}.log"
    done
done
#!/bin/bash

## Remove SNPs located in low complexity and simple repeat regions
echo "Removing SNPs in low complexity and simple repeat regions"
# Step 1: Convert rmsk file to BED format
convert2bed -i rmsk <  /home/roo/data/filtersnp/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_rm.out >  /home/roo/data/filtersnp/out.bed

# Step 2: Extract lines containing "Simple_repeat" and "Low_complexity" into separate files
grep "Simple_repeat"  /home/roo/data/filtersnp/out.bed >  /home/roo/data/filtersnp/repeat.bed
grep "Low_complexity"  /home/roo/data/filtersnp/out.bed >  /home/roo/data/filtersnp/low.bed

# Combine low.bed into repeat.bed
cat  /home/roo/data/filtersnp/low.bed >>  /home/roo/data/filtersnp/repeat.bed

# Step 3: Extract only the first three columns and save to lowrepeat.bed
awk '{print $1"\t"$2"\t"$3}'  /home/roo/data/filtersnp/repeat.bed >  /home/roo/data/filtersnp/lowrepeat.bed

# Step 4: Run vcftools command
vcftools --vcf ~/data/raw_vcf/ALL.snp.filter.vcf --exclude-bed  /home/roo/data/filtersnp/lowrepeat.bed --out rmlowrepeat --recode

## Remove sex scaffolds
echo "Removing sex scaffolds"
#1. Detect sex scaffolds using the SATC method with the satcFunc.R (https://github.com/popgenDK/SATC/tree/main) function, generating sexscaffolds.txt

#2. Remove SNPs located on sex scaffolds from the VCF file:
awk 'NR==FNR {scaffold[$1]; next} !($1 in scaffold)'  /home/roo/data/filtersnp/sexscaffolds.txt  /home/roo/data/filtersnp/rmlowrepeat.recode.vcf >  /home/roo/data/filtersnp/rmlowrepeat.sex.vcf

## Filter SNPs based on sequencing depth, missing rate, and minor allele frequency
echo "Filtering SNPs based on sequencing depth, missing rate, and minor allele frequency"
# 1. Calculate sequencing depth for each sample
for sample in BS001 BS002 BS003 BS004 BS005 BS006 BS007 BS008 BT001 BT002 BT003 BT004 DT001 DT002 DT003 DT004 DT005 DT006 DT007 DT008 DT009 DT010 DT011 DT012 DT013 DT014 DT015 DT017 DT018 DT019 DT020 DT021 DT022 DT023 DT024 DT025 DT027 DT028 DT029 DT030 DT031 FA006 FA021 FA028 FA029 FA030 FA032 FA033 FA036 HK004 HK005 LS001 QH001 QH002 QH003 QH004 QH005 QH006 QH007 QH008 QH009 QH010 QH011 QH012 QH013 QH014 QH015 QH016 QH017 QH018 QH019 QH020 QH021 SY003 SY004 SY005 SY006 SY007 SY008 WC033 WC035 WC038 WC039 WC040 WC041 WC042 FA007 FA012 FA013 FA015 FA016 FA019 FA020 FA024 FA026 FA035 FA039 HK003 SY001 WC008 WC009 WC011 WC013 WC032 WC034 WC036 WC037 WC043 WC044 WC045
do
    /home/roo/data/filtersnp/PanDepth-2.25-Linux-x86_64/pandepth -i /home/data/qj/BAM/bam/$sample.rmdup.bam -o /home/roo/data/filtersnp/sample/$sample
done

# 2. Filter biallelic SNPs with depth between 10 and 40
vcftools --maf 0.05 --max-missing 0.8 --min-alleles 2 --max-alleles 2 --minDP 10 --maxDP 40 --vcf  /home/roo/data/filtersnp/rmlowrepeat_sex.vcf --recode --out  /home/roo/data/filtersnp/vcfrmlowrepeat.sex.depth.mismafali

## Remove closely related individuals
echo "Removing closely related individuals"
# Step 1: Convert VCF file to PLINK format
plink2 --allow-extra-chr --vcf /home/roo/data/filtersnp//home/roo/data/filtersnp/vcfrmlowrepeat.sex.depth.mismafali.recode.vcf --make-bed --out kinship

# Step 2: Calculate relatedness and generate PLINK files
plink2 --allow-extra-chr --bfile kinship --king-cutoff 0.177 --make-bed --out kinship

# Step 3: Extract IDs of individuals to remove from kinship.king.cutoff.out.id file (excluding header)
awk 'NR > 1 {print $2}' kinship.king.cutoff.out.id > samples_to_remove.txt

# Step 4: Remove closely related individuals from the original VCF file
vcftools --vcf /home/roo/data/filtersnp/vcfrmlowrepeat.sex.depth.mismafali.recode.vcf --remove samples_to_remove.txt --recode --out /home/roo/data/filtersnp/datset1

## LD Filtering
echo "LD filtering"
# Step 1: Convert VCF file to PLINK format
plink --allow-extra-chr --vcf /home/roo/data/filtersnp/datset1.recode.vcf --make-bed --out filtered

# Step 2: Perform LD filtering
plink --allow-extra-chr --bfile filtered --indep-pairwise 50 5 0.2 --out dataset2

# Step 3: Apply LD filtering
plink --allow-extra-chr --bfile filtered --extract dataset.prune.in --make-bed --out dataset2

# Step 4: Convert PLINK format back to VCF format
plink --allow-extra-chr --bfile dataset2 --recode vcf --out dataset2
#!/bin/bash

plink --allow-extra-chr --vcf /home/roo/data/anlaysis_revise/dataset1/dataset1.recode.vcf --make-bed --out dataset1
plink --allow-extra-chr --bfile /home/roo/data/anlaysis_revise/dataset1/dataset1 --out dataset1.pc10 --pca 10 --threads 20
#!/bin/bash
# 1. Format conversion
python vcf2phylip.py --input /home/roo/data/anlaysis_revise/dataset1/dataset1.recode.vcf

# 2. Build IQ-TREE
iqtree2 -s /home/roo/data/anlaysis_revise/dataset1/ML/dataset1.recode.min4.phy -st DNA -T 2 -m GTR+ASC -redo -B 1000 -bnni --prefix iqtree
#!/bin/bash

# Set variables
REFERENCE="/home/roo/data/pygenome/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.fna"  # Reference genome
INPUT_DIR="/home/data/qj/test"                   # Directory with input FASTQ files
OUTPUT_DIR="/home/roo/data/raw_vcf"                 # Output directory
SAMPLES=("BS001" "BS002" "BS003" "BS004" "BS005" "BS006" "BS007" 
         "BS008" "BT001" "BT002" "BT003" "BT004" "DT001" "DT002" 
		 "DT003" "DT004" "DT005" "DT006" "DT007" "DT008" "DT009" 
		 "DT010" "DT011" "DT012" "DT013" "DT014" "DT015" "DT017" 
		 "DT018" "DT019" "DT020" "DT021" "DT022" "DT023" "DT024" 
		 "DT025" "DT027" "DT028" "DT029" "DT030" "DT031" "FA012" 
		 "FA013" "FA015" "FA016" "FA019" "FA020" "FA021" "FA024" 
		 "FA026" "FA028" "FA029" "FA030" "FA032" "FA033" "FA035" 
		 "FA036" "FA039" "HK004" "HK005" "LS001" "QH001" "QH002" 
		 "QH003" "QH004" "QH005" "QH006" "QH007" "QH008" "QH009" 
		 "QH010" "QH011" "QH012" "QH013" "QH014" "QH015" "QH016" 
		 "QH017" "QH018" "QH019" "QH020" "QH021" "SY003" "SY004" 
		 "SY005" "SY006" "SY007" "SY008" "WC008" "WC009" "WC011" 
		 "WC013" "FA006" "FA007" "WC032" "WC033" "WC034" "WC035" 
		 "WC036" "WC037" "WC038" "WC039" "WC040" "WC041" "WC042" 
		 "WC043" "WC044" "WC045" "ZO001" "ZO002") # Sample names list

# Create output directory
mkdir -p ${OUTPUT_DIR}

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}"

    INPUT_FASTQ1="${INPUT_DIR}/${SAMPLE}_1.fq.gz"
    INPUT_FASTQ2="${INPUT_DIR}/${SAMPLE}_2.fq.gz"
	clean_fastq1="${INPUT_DIR}/${SAMPLE}_1_clean.fq.gz"
	clean_fastq2="${INPUT_DIR}/${SAMPLE}_2_clean.fq.gz"
    BAM_FILE="${OUTPUT_DIR}/${SAMPLE}.bam"
    SORTED_BAM_FILE="${OUTPUT_DIR}/${SAMPLE}_sorted.bam"
    RECAL_BAM_FILE="${OUTPUT_DIR}/${SAMPLE}_recal.bam"
	GVCF_FILE="${OUTPUT_DIR}/${SAMPLE}.g.vcf"
    RECAL_DATA_FILE="${OUTPUT_DIR}/${SAMPLE}_recal_data.table"
	
	
	# 1. Quality assessment
	echo "1.Quality assessment"
	mkdir fastqc_result #Directory for quality assessment results
    fastqc ${INPUT_FASTQ1}
	fastqc ${INPUT_FASTQ2}
    mv *.zip *.html fastqc_result/
	
	# 2. Sequencing data filtering
	echo "2.Sequencing data filtering"
    mkdir fastp_out #Directory for filtered results
    fastp -i ${INPUT_FASTQ1} ${INPUT_FASTQ2} -o ${clean_fastq1}  -O ${clean_fastq2} -z 4 -q 20 -u 30 -n 10 -f 15 -t 15 -F15 -T 15
	
    # 3. Alignment
    echo "Step 3: Aligning reads to the reference genome"
    bwa index ${REFERENCE}
    bwa mem ${REFERENCE} ${INPUT_FASTQ1} ${INPUT_FASTQ2} | samtools view -bS - > ${BAM_FILE}

    # 4. Sorting
    echo "Step 4: Sorting BAM file"
    samtools sort ${BAM_FILE} -o ${SORTED_BAM_FILE}

    # 5.  Indexing
    echo "Step 5: Indexing BAM file"
    samtools index ${SORTED_BAM_FILE}
	
	# 6. Statistics for alignment and sequencing depth
	echo "6.Statistics for alignment and sequencing depth"
    samtools depth ${SORTED_BAM_FILE} >${SORTED_BAM_FILE}_depth.txt #将测序深度保存为depth.txt
    samtools flagstat ${SORTED_BAM_FILE} >${SORTED_BAM_FILE}_mapping.txt #将比对结果保存为T_result.txt

    # 7. Recalibration (using GATK)
    echo "Step 7: Recalibration (optional)"
    gatk MarkDuplicates -I ${SORTED_BAM_FILE} -O ${RECAL_BAM_FILE} -M ${RECAL_DATA_FILE}
    
	# 8. Generate a .dict file for the reference genome
    gatk CreateSequenceDictionary -R ${REFERENCE} -O ${REFERENCE}.dict
    
	# 9. Generate .g.vcf file for each sample
	Gatk HaplotypeCaller –R ${REFERENCE} –emit-ref-confidence GVCF –I ${RECAL_DATA_FILE} –O ${GVCF_FILE}
	
	# 10.Merge g.vcf files
    find ${OUTPUT_DIR} -name "*.g.vcf" > input.list
    gatk CombineGVCFs  -R ${REFERENCE} --variant input.list -O combine.g.vcf

    # 11.Multi-sample variant calling
	echo "Step 11: Calling variants"
    gatk GenotypeGVCFs  -R ${REFERENCE}  -V combined.g.vcf  -O raw.vcf

    # 12. Filter variants
    echo "Step 12: Filtering variants"
	gatk VariantFiltration -V raw.vcf  --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || QUAL < 30.0" 
	--filter-name "PASS" -O raw.filter.vcf
	
	#13.Select SNPs
    gatk SelectVariants  -select-type SNP -V raw.filter.vcf -O ALL.snp.filter.vcf
    
done

echo "All samples have been processed."
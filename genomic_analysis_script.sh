#!/bin/bash
set -e

SRA=SRR1972739
REF_ID=AF086833.2
RESULTS_FOLDER=results
RAW_DIR=${RESULTS_FOLDER}/raw
ALIGNED_DIR=${RESULTS_FOLDER}/aligned
VARIANT_DIR=${RESULTS_FOLDER}/variants
ANNOTATED_DIR=${RESULTS_FOLDER}/annotated
QC_DIR=${RESULTS_FOLDER}/qc
SNPEFF_DIR=${RESULTS_FOLDER}/snpEff
SNPEFF_DATA_DIR=${SNPEFF_DIR}/data/reference_db

mkdir -p $RAW_DIR $ALIGNED_DIR $VARIANT_DIR $ANNOTATED_DIR $QC_DIR $SNPEFF_DATA_DIR

echo Downloading reference genome...
efetch -db nucleotide -id $REF_ID -format fasta > $RAW_DIR/reference.fasta
echo Downloaded reference genome!

echo Downloading sequencing data...
prefetch $SRA -O $RAW_DIR
fastq-dump -X 10000 $RAW_DIR/${SRA}/${SRA}.sra -O $RAW_DIR
echo Downloaded sequencing data!

if [ ! -s $RAW_DIR/reference.fasta ]; then
  echo Error: Reference genome file is missing or empty. >&2
  exit 1
fi

if [ ! -s $RAW_DIR/${SRA}.fastq ]; then
  echo Error: FASTQ file is missing or empty. >&2
  exit 1
fi

echo Running FastQC on raw reads...
fastqc -o $QC_DIR $RAW_DIR/${SRA}.fastq

echo Indexing reference genome with samtools...
samtools faidx $RAW_DIR/reference.fasta

echo Building BWA index...
bwa index $RAW_DIR/reference.fasta

echo Creating FASTA dictionary using GATK...
gatk CreateSequenceDictionary -R $RAW_DIR/reference.fasta -O $RAW_DIR/reference.dict

echo Aligning reads with read groups...
bwa mem -R '@RG\tID:1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:sample1' $RAW_DIR/reference.fasta $RAW_DIR/${SRA}.fastq > $ALIGNED_DIR/aligned.sam
echo Aligned reads!

echo Converting SAM to sorted BAM...
samtools view -b $ALIGNED_DIR/aligned.sam | samtools sort -o $ALIGNED_DIR/aligned.sorted.bam

echo Validating BAM file...
gatk ValidateSamFile -I $ALIGNED_DIR/aligned.sorted.bam -MODE SUMMARY

echo Marking duplicates...
gatk MarkDuplicates -I $ALIGNED_DIR/aligned.sorted.bam -O $ALIGNED_DIR/dedup.bam -M $ALIGNED_DIR/dup_metrics.txt

echo Indexing deduplicated BAM file...
samtools index $ALIGNED_DIR/dedup.bam

echo Calling variants...
gatk HaplotypeCaller -R $RAW_DIR/reference.fasta -I $ALIGNED_DIR/dedup.bam -O $VARIANT_DIR/raw_variants.vcf
echo Called variants!

echo Filtering variants...
gatk VariantFiltration -R $RAW_DIR/reference.fasta -V $VARIANT_DIR/raw_variants.vcf -O $VARIANT_DIR/filtered_variants.vcf --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER

echo Downloading reference GenBank file for snpEff...
efetch -db nucleotide -id $REF_ID -format genbank > $SNPEFF_DATA_DIR/genes.gbk
echo Downloaded GenBank file for snpEff!

echo Creating custom snpEff configuration file...
cat <<EOF > $SNPEFF_DIR/snpEff.config
# Custom snpEff config for reference_db
reference_db.genome : reference_db
reference_db.fa : $(readlink -f $RAW_DIR/reference.fasta)
reference_db.genbank : $(readlink -f $SNPEFF_DATA_DIR/genes.gbk)
EOF

echo Building snpEff database...
snpEff build -c $SNPEFF_DIR/snpEff.config -genbank -v -noCheckProtein reference_db
echo Built snpEff database!

echo Exporting snpEff database...
snpEff dump -c $SNPEFF_DIR/snpEff.config reference_db > $SNPEFF_DIR/snpEff_reference_db.txt 
echo Exported snpEff database!

echo Annotating variants with snpEff...
snpEff -c $SNPEFF_DIR/snpEff.config -stats $SNPEFF_DIR/snpEff.html reference_db $VARIANT_DIR/filtered_variants.vcf > $ANNOTATED_DIR/annotated_variants.vcf
echo Annotated variants with snpEff!

echo Pipeline completed successfully! Check the folders in $RESULTS_FOLDER for output files.
tree
SRA="SRR1972739"
REF_ID="AF086833.2"
RESULTS_FOLDER="results"
RAW_DIR=f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR=f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR=f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR=f"{RESULTS_FOLDER}/annotated"
QC_DIR=f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR=f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR=f"{SNPEFF_DIR}/data/reference_db"
SNAKEMAKE_DIR=f"{RESULTS_FOLDER}/snakemake"

rule all:
    input:
        f"{SNAKEMAKE_DIR}/dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}.fastq",
        f"{QC_DIR}/{SRA}_fastqc.html",
        f"{QC_DIR}/{SRA}_fastqc.zip",
        f"{RAW_DIR}/reference.dict",
        f"{ALIGNED_DIR}/aligned.sam",
        f"{ALIGNED_DIR}/aligned.sorted.bam",
        f"{ALIGNED_DIR}/dedup.bam",
        f"{ALIGNED_DIR}/dup_metrics.txt",
        f"{VARIANT_DIR}/raw_variants.vcf",
        f"{VARIANT_DIR}/filtered_variants.vcf",
        f"{SNPEFF_DATA_DIR}/genes.gbk",
        f"{SNPEFF_DIR}/snpEff.config",
        f"{ANNOTATED_DIR}/annotated_variants.vcf",
        f"{SNPEFF_DIR}/snpEff.html"
        

rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/dirs_created"
    shell:
        """
        echo "Making Pipeline files..."
        mkdir -p {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR}
        touch {output.marker}
        """

rule download_fasta:
    input:
        f"{SNAKEMAKE_DIR}/dirs_created"    
    output:
        f"{RAW_DIR}/reference.fasta" 
    shell:
        """
        echo Downloading reference genome...
        efetch -db nucleotide -id {REF_ID} -format fasta > {RAW_DIR}/reference.fasta
        echo Downloaded reference genome!
        """

rule download_fastq:
    input:
        f"{SNAKEMAKE_DIR}/dirs_created"
        
    output:
        f"{RAW_DIR}/{SRA}.fastq"
    
    shell:
        """
        echo Downloading sequencing data...
        prefetch {SRA} -O {RAW_DIR}
        fastq-dump -X 10000 {RAW_DIR}/{SRA}/{SRA}.sra -O {RAW_DIR}
        echo Downloaded sequencing data!
        """

rule QC_check:
    input:
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}.fastq"
        
    output:
        f"{QC_DIR}/{SRA}_fastqc.html",
        f"{QC_DIR}/{SRA}_fastqc.zip"

    shell:
        """
        if [ ! -s {RAW_DIR}/reference.fasta ]; then
            echo Error: Reference genome file is missing or empty. >&2
            exit 1
        fi

        if [ ! -s {RAW_DIR}/{SRA}.fastq ]; then
            echo Error: FASTQ file is missing or empty. >&2
            exit 1
        fi

        echo Running FastQC on raw reads...
        fastqc -o {QC_DIR} {RAW_DIR}/{SRA}.fastq
        """

rule Create_SAM:
    input:
        f"{SNAKEMAKE_DIR}/dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}.fastq"

    output:
        f"{RAW_DIR}/reference.dict",
        f"{ALIGNED_DIR}/aligned.sam"
        
    shell:
        """
        echo Indexing reference genome with samtools...
        samtools faidx {RAW_DIR}/reference.fasta

        echo Building BWA index...
        bwa index {RAW_DIR}/reference.fasta

        echo Creating FASTA dictionary using GATK...
        gatk CreateSequenceDictionary -R {RAW_DIR}/reference.fasta -O {RAW_DIR}/reference.dict

        echo Aligning reads with read groups...
        bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:sample1' {RAW_DIR}/reference.fasta {RAW_DIR}/{SRA}.fastq > {ALIGNED_DIR}/aligned.sam
        echo Aligned reads!
        """

rule Create_BAM:
    input:
        f"{ALIGNED_DIR}/aligned.sam"

    output:
        f"{ALIGNED_DIR}/aligned.sorted.bam",
        f"{ALIGNED_DIR}/dedup.bam",
        f"{ALIGNED_DIR}/dedup.bam.bai",
        f"{ALIGNED_DIR}/dup_metrics.txt"

    shell:
        """
        echo Converting SAM to sorted BAM...
        samtools view -b {input}| samtools sort -o {output[0]}

        echo Validating BAM file...
        gatk ValidateSamFile -I {output[0]} -MODE SUMMARY

        echo Marking duplicates...
        gatk MarkDuplicates -I {output[0]} -O {output[1]} -M {output[3]}

        echo Indexing deduplicated BAM file...
        samtools index {output[1]}
        """

rule Variant_Calling:
    input:
        f"{RAW_DIR}/reference.fasta",
        f"{ALIGNED_DIR}/dedup.bam"
    output:
        f"{VARIANT_DIR}/raw_variants.vcf",
        f"{VARIANT_DIR}/filtered_variants.vcf"
    shell:
        """
        echo Calling variants...
        gatk HaplotypeCaller -R {RAW_DIR}/reference.fasta -I {ALIGNED_DIR}/dedup.bam -O {VARIANT_DIR}/raw_variants.vcf
        echo Called variants!

        echo Filtering variants...
        gatk VariantFiltration -R {RAW_DIR}/reference.fasta -V {VARIANT_DIR}/raw_variants.vcf -O {VARIANT_DIR}/filtered_variants.vcf --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER
        echo Variants filtered!
        """

rule snpEff_database:
    input:
        f"{RAW_DIR}/reference.fasta"
    output:
        f"{SNPEFF_DATA_DIR}/genes.gbk",
        f"{SNPEFF_DIR}/snpEff.config",
        f"{SNPEFF_DIR}/snpeff_build.done",
        f"{SNPEFF_DIR}/snpEff_reference_db.txt",

    shell:
        """
        echo Downloading reference GenBank file for snpEff...
        efetch -db nucleotide -id {REF_ID} -format genbank > {SNPEFF_DATA_DIR}/genes.gbk
        echo Downloaded GenBank file for snpEff!

        echo Creating custom snpEff configuration file...
        cat <<EOF > {SNPEFF_DIR}/snpEff.config
# Custom snpEff config for reference_db
reference_db.genome : reference_db
reference_db.fa : $(readlink -f {RAW_DIR}/reference.fasta)
reference_db.genbank : $(readlink -f {SNPEFF_DATA_DIR}/genes.gbk)
EOF

        echo Building snpEff database...
        snpEff build -c {SNPEFF_DIR}/snpEff.config -genbank -v -noCheckProtein reference_db
        touch {SNPEFF_DIR}/snpeff_build.done
        echo Built snpEff database!

        echo Exporting snpEff database...
        snpEff dump -c {SNPEFF_DIR}/snpEff.config reference_db > {SNPEFF_DIR}/snpEff_reference_db.txt
        echo Exported snpEff database!
        """

rule snpEff_Annotation:
    input:
        f"{SNPEFF_DIR}/snpEff.config",
        f"{VARIANT_DIR}/filtered_variants.vcf"
    output:
        f"{ANNOTATED_DIR}/annotated_variants.vcf",
        f"{SNPEFF_DIR}/snpEff.html"
    shell:
        """
        echo Annotating variants with snpEff...
        snpEff -c {SNPEFF_DIR}/snpEff.config -stats {SNPEFF_DIR}/snpEff.html reference_db {VARIANT_DIR}/filtered_variants.vcf > {ANNOTATED_DIR}/annotated_variants.vcf
        echo Annotated variants with snpEff!
        echo Pipeline completed successfully! Check the folders in {SNPEFF_DIR} for output files.
        """

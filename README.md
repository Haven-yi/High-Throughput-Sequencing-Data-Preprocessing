# High-Throughput Sequencing Data Preprocessing
This tutorial provides a step-by-step guide for preprocessing high-throughput sequencing data in an HPC environment using Slurm array jobs. The workflow utilizes Fastp for quality trimming, followed by alignment with BWA-MEM2, sorting with Samtools, and duplicate removal with Picard, ensuring efficient and accurate preprocessing for variant calling.

Required Tools and Files:  
Fastp (for quality trimming)  
BWA-MEM2 (for read alignment)  
Samtools (for BAM file sorting and statistics)  
Picard (for duplicate removal)  

Required Files:  
Reference Genome File: genome.fasta (FASTA format containing the reference genome sequences)  
Sequencing Data Files: Paired-end FASTQ files, such as S1_1.fq.gz and S1_2.fq.gz.

# Step 1: Quality Trimming
```
#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**
#SBATCH --mem=**
#SBATCH --time=**
#SBATCH --array=1-**
#SBATCH --output=logs/fastp_%A_%a.out
#SBATCH --error=logs/fastp_%A_%a.err

module load StdEnv/2023
module load fastp/0.24.0
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample.txt)
R1="${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE}_R2.fastq.gz"
CLEAN_READS_DIR="/home/fasta-cleanreads"
REPORT_DIR="/home/fasta-report"
SAMPLE_NAME=$(basename "$SAMPLE")
CLEAN_R1="${CLEAN_READS_DIR}/${SAMPLE_NAME}_clean_read1.fq.gz"
CLEAN_R2="${CLEAN_READS_DIR}/${SAMPLE_NAME}_clean_read2.fq.gz"
REPORT_HTML="${REPORT_DIR}/${SAMPLE_NAME}_report.html"
REPORT_JSON="${REPORT_DIR}/${SAMPLE_NAME}_report.json"

fastp -i "$R1" -I "$R2" -o "$CLEAN_R1" -O "$CLEAN_R2" \
    --trim_front1 10 --trim_front2 10 \
    --trim_tail1 10 --trim_tail2 10 \
    --cut_front --cut_tail \
    --qualified_quality_phred 20 \
    --detect_adapter_for_pe \
    --thread ** \
    --html "$REPORT_HTML" \
    --json "$REPORT_JSON"
```

# Step 2: Build Genome Index
```
samtools faidx genome.fasta
bwa-mem2 index genome.fasta
samtools dict genome.fasta -o genome.dict
```
# Step 3: Align Reads to the Reference Genome
Once the genome index is built, the next step is to align the sequencing reads to the reference genome using BWA-MEM2. After alignment, we will use Samtools to convert the resulting SAM file to BAM format and sort it.
```
#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**
#SBATCH --mem=**
#SBATCH --time=**
#SBATCH --array=1-**
#SBATCH --output=bwalogs/bwa_%A_%a.out
#SBATCH --error=bwalogs/bwa_%A_%a.err

module load StdEnv/2023
module load bwa-mem2/2.2.1
module load gcc/12.3
module load samtools/1.20

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample.txt)
SAMPLE_NAME=$(basename "$SAMPLE")
R1="${SAMPLE}_clean_read1.fq.gz"
R2="${SAMPLE}_clean_read2.fq.gz"
REF="/home/genome.fasta"
OUTPUT_BAM="/home/bam/${SAMPLE_NAME}.sort.bam"
LOG="/home/bam/${SAMPLE_NAME}.bwa.log"

bwa-mem2 mem -t 8 -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:illumina" \
    "$REF" "$R1" "$R2" 2> "$LOG" | \
samtools sort -@ ** -m ** -o "$OUTPUT_BAM"
```

# Step 4: Remove Duplicates
Duplicate reads may arise due to PCR amplification and should be removed to avoid bias in variant calling. 
```
#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**
#SBATCH --mem=**
#SBATCH --time=**
#SBATCH --array=1-**
#SBATCH --output=markdup-logs/markdup_%A_%a.out
#SBATCH --error=markdup-logs/markdup_%A_%a.err

NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample.txt)
SAMPLE="/home/bam/${NAME}.sort.bam"
MARKDUPBAM="/home/markdup-bam/$NAME.markdup.sort.bam"
MARKDUP="/home/markdup/$NAME.metrics"

module load picard
module load StdEnv/2023
module load samtools

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT=$SAMPLE \
    OUTPUT=/dev/stdout \
    METRICS_FILE=$MARKDUP \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT | \
samtools sort -@ ** -o $MARKDUPBAM
```


# Step 5: Quality Control and Statistics
After removing duplicates, it's important to check the quality of the BAM files. This includes alignment statistics and coverage information.
```
samtools flagstat S1.sorted.bam > S1.sorted.bam.flagstat
samtools coverage S1.sorted.bam > S1.sorted.bam.coverage
```
These steps prepare your data for downstream variant calling and ensure that only high-quality data is used for further analysis.

# _M. tuberculosis_ Whole Genome Sequence Processing Workflow

## 1. Schematic Representation for the Workflow
![Image of Workflow](https://github.com/cwchang-nereus/Mtb_WGS_Processing/blob/main/WGS_workflow.png)

## 2. Environment Settings for Sequence Data Processing
### 2.1 Conda session
#### 2.1.1 Create a conda environment for reads quality control
```
conda create -n mtb_wgs_qc_env
```

### 2.2 Install WGS-analysing tools 
#### 2.2.1 Quality control
```
conda install \
  --freeze-installed \
  -n mtb_wgs_qc_env \
  -c bioconda \
  fastqc trimmomatic fastp multiqc qualimap
```

### 2.3 Create a conda environment for reads mapping and variant calling
```
conda create -n mtb_wgs_map_vc_env
```

### 2.4 Reads and variants processing and filtering tools
#### 2.4.1 install ncurses rom conda-forge before samtools
```
conda install \
  -c conda-forge \
  ncurses
```
```
conda install \
  --freeze-installed \
  -n mtb_wgs_map_vc_env \
  -c bioconda \
  delly bwa freebayes bowtie2 samtools bedtools bcftools vcftools matplotlib
```

### 2.5 Create a conda environment for reads de novo assembly and variant calling
```
conda create -n mtb_wgs_denovo_vc_env
```

#### 2.5.1 de nova assembly tools
```
conda install \
  --freeze-installed \
  -n mtb_wgs_denovo_vc_env \
  -c bioconda \
  velvet spades samtools bedtools bcftools vcftools
```

### 2.6 Create a conda environment for data visualization
```
conda create -n mtb_wgs_vis_env
```

#### 2.6.1 Genomic data visualization tools
```
conda install \
  --freeze-installed \
  -n mtb_wgs_vis_env \
  -c bioconda \
  igv igvtools
```

### 2.7 Create a conda environment for data retrieving
```
conda create -n data_retrv_env
```

#### 2.7.1 Databases connection
```
conda install \
  --freeze-installed \
  -n data_retrv_env \
  -c bioconda \
  ncbi-genome-download
```

### 2.8 Create a conda environment for vcf annotation
```
conda create -n mtb_wgs_anno_env
```

#### 2.8.1 annotation and analysis tools
```
conda install \
  --freeze-installed \
  -n mtb_wgs_anno_env \
  -c bioconda \
  snpeff snpsift bcftools snp-dists seqmagick
```

### 2.9 Local session
#### 2.9.1 Brew session
```
brew install bzip2 pbzip2 gzip
brew install python python@2
brew install picard-tools
brew install git-lfs
```

#### 2.9.2 GATK
##### 2.9.2.1 GATK 4.0
```
TOOLS_PATH=[PATH to WHERE TOOL FILES LOCATED]
GATK_PATH=$TOOLS_PATH/GATK

cd $GATK_PATH
mkdir GATK4

git clone https://github.com/broadinstitute/gatk.git
mv gatk ./GATK4/

cd GATK4
git lfs pull
git lfs install
./gradlew bundle
```

#### 2.9.2.2 GATK 3.8-1-0
```
mkdir -p $GATK_PATH/GATK3.8-1-0
wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef" \
  -O $GATK_PATH/GATK3.8-1-0/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

bunzip2 GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
tar -xvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar
```

#### 2.9.2.3 Set alias for GATK
```
alias gatk="java -jar -Xmx64g $(find $GATK_PATH/GATK4 -name '*local.jar')"
alias gatk3="java -jar -Xmx64g $(find $GATK_PATH/GATK3.8-1-0 -name '*.jar')"
```

#### 2.9.3 Choose specific version of Java in Ubuntu
```
sudo update-alternatives --config java
```

#### 2.9.4Alfred
```
docker pull trausch/alfred
```


## 3. Data Preparation
### 3.1 Compress raw fastq files to .gz files
parallel compressing
```
WORKING_PATH=[PATH to A WORKING DIREACTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

FASTQF_PATH=$WORKING_PATH/01_SeqData/rawReads
FASTQ_FILE=$(ls $FASTQF_PATH)
FASTQGZ_PATH=$WORKING_PATH/01_SeqData/rawReads.gz

cd $FASTQF_PATH

pigz -6 $(echo $FASTQ_FILE) -kv -p 4
mv ./*.fastq.gz $FASTQGZ_PATH 
```

### 3.2 Reference Sequences Downloading
Download the reference sequences of _M. tuberculosis_
```
conda activate data_retrv_env

REF_GENOME_PATH=$REF_PATH/11_Reference_Genome/

mkdir -p $REF_GENOME_PATH

ncbi-genome-download \
  --format "all" \
  --assembly-level complete \
  --genus 'Mycobacterium tuberculosis' \
  --refseq-category reference \
  --output-folder $REF_GENOME_PATH \
  bacteria

# H37Rv
cd $REF_GENOME_PATH/refseq/bacteria/GCF_000195955.2

for GZ in feature_table.txt.gz cds_from_genomic.fna.gz genomic.fna.gz genomic.gbff.gz genomic.gff.gz
  do
    gunzip -cv \
    $REF_GENOME_PATH/refseq/bacteria/GCF_000195955.2/GCF_000195955.2_ASM19595v2_$GZ \
    > $REF_GENOME_PATH/refseq/bacteria/GCF_000195955.2/NC_000962.3_Mtb_H37Rv.${GZ/.gz/}
  done

# bovis
cd $REF_GENOME_PATH/refseq/bacteria/GCF_000195835.2

for GZ in feature_table.txt.gz cds_from_genomic.fna.gz genomic.fna.gz genomic.gbff.gz genomic.gff.gz
  do
    gunzip -cv \
    $REF_GENOME_PATH/refseq/bacteria/GCF_000195835.2/GCF_000195835.2_ASM19583v2_$GZ \
    > $REF_GENOME_PATH/refseq/bacteria/GCF_000195835.2/NC_002945.4_Mbovis.${GZ/.gz/}
  done
```
### 3.3 Raw reads quality control
#### 3.3.1 Generate raw read quality reports
FASTQ report
```
conda activate mtb_wgs_qc_env

WORKING_PATH=[PATH to A WORKING DIREACTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

RAW_FASTQGZ_PATH=$WORKING_PATH/01_SeqData/rawReads.gz
RAW_FASTQ_OUTPUT_PATH=$WORKING_PATH/02_Quality_Control/raw_read_fastq_report

mkdir -p $RAW_FASTQ_OUTPUT_PATH
cd $RAW_FASTQGZ_PATH

for FASTQGZ_FILE in $(ls $RAW_FASTQGZ_PATH) 
  do 
    fastqc $FASTQGZ_FILE \
    --outdir $RAW_FASTQ_OUTPUT_PATH \
    --threads 16
  done
```

#### 3.3.2 Raw reads quality trimming
by FASTP
```
RAW_READ_QC_OUTPUT_PATH=$WORKING_PATH/02_Quality_Control/raw_read_fastp_qc_report

mkdir -p $RAW_READ_QC_OUTPUT_PATH
cd $RAW_FASTQGZ_PATH

for READ1 in $(ls | grep R1)
  do
    READ2=${READ1/R1/R2}
    REPORT_ID=$(echo $READ1 | cut -d"_" -f 1)
    
    fastp \
      --in1 $READ1 \
      --out1 $RAW_READ_QC_OUTPUT_PATH/${REPORT_ID}_trimmed_filtered_R1.fastq.gz \
      --in2 $READ2 \
      --out2 $RAW_READ_QC_OUTPUT_PATH/${REPORT_ID}_trimmed_filtered_R2.fastq.gz \
      --overrepresentation_analysis  `# overrepresented sequence analysis`  \
      --correction                   `# base correction in overlapped regions (only for PE data)`  \
      --overlap_len_require 30       `# the minimum length to detect overlapped region of PE reads.` \
      --detect_adapter_for_pe        `# detection for adapter on both read1 and read2`  \
      --qualified_quality_phred 30   `# set the quality value that a base is qualified. Q>=Q30 is qualified` \
      --length_required 30           `# reads shorter than 30 will be discarded` \
      --cut_tail                     `# drop the bases in the window if its mean quality < threshold` \
      --cut_tail_window_size 4       `# the window size option of cut_tail` \
      --cut_tail_mean_quality 20     `# the mean quality requirement option for cut_tail` \
      --html $RAW_READ_QC_OUTPUT_PATH/${REPORT_ID}_trimmed_filtered_fastp.html \
      --json $RAW_READ_QC_OUTPUT_PATH/${REPORT_ID}_trimmed_filtered_fastp.json \
      --thread 16
  done
```

#### 3.3.3 FASTQ report generation
```
QCED_FASTQ_OUTPUT_PATH=$WORKING_PATH/02_Quality_Control/qced_read_fastq_report
RAW_READ_QC_OUTPUT_PATH=$WORKING_PATH/02_Quality_Control/raw_read_fastp_qc_report

mkdir -p $QCED_FASTQ_OUTPUT_PATH
cd $RAW_READ_QC_OUTPUT_PATH

for FASTQGZ_FILE in $(ls $RAW_READ_QC_OUTPUT_PATH | grep fastq.gz) 
  do 
    fastqc $FASTQGZ_FILE \
    --outdir $QCED_FASTQ_OUTPUT_PATH \
    --threads 16
  done
```

# 4.Reads Mapping
## 4.1 Alingment by BWA-MEM algorithm
### 4.1.1 Index the reference genome and create the dictionary
```
conda activate mtb_wgs_map_vc_env

WORKING_PATH=[PATH to A WORKING DIREACTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

REF_GENOME_PATH=$REF_PATH/11_Reference_Genome/

samtools faidx $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna")

gatk CreateSequenceDictionary \
  --REFERENCE $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --OUTPUT $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna" | sed 's/.fna/.dict/')
```

### 4.1.2 Mapping sample reads to a reference genome by bwa-mem
```
RAW_READ_QC_OUTPUT_PATH=$WORKING_PATH/02_Quality_Control/raw_read_fastp_qc_report
SAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/01_SAM

mkdir -p $SAM_FILE_PATH
cd $RAW_READ_QC_OUTPUT_PATH

for READ1 in $(find $RAW_READ_QC_OUTPUT_PATH -name "*R1*.fastq.gz" -exec basename "{}" \;)
  do
    READ2=${READ1/R1/R2}
    SAMFILE_ID=$(echo $READ1 | cut -d"_" -f 1)
    
    bwa mem \
      -t 16 \
      $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      $READ1 $READ2 \
      > $SAM_FILE_PATH/${SAMFILE_ID}.sam
  done
```

### 4.1.3 Clean up read pairing information and flags left by BWA, sort into coordinate order and index BAM files
```
SAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/01_SAM
BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/02_BAM

mkdir -p $BAM_FILE_PATH

cd $SAM_FILE_PATH

for SAMFILE in $(find $SAM_FILE_PATH -name "*.sam" -exec basename "{}" \;)
  do
    printf "\n>>> initiating fixmate..."
    
    samtools fixmate \
      -@ 16 \
      $SAMFILE \
      $BAM_FILE_PATH/${SAMFILE/.sam/_fixmate.bam} \
      --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      --output-fmt BAM
    
    printf "\n>>> fixmate done: ${SAMFILE/.sam/_fixmate.bam} generated"
    printf "\n>>> initiating sort..."
    
    samtools sort \
      -@ 16 \
      $BAM_FILE_PATH/${SAMFILE/.sam/_fixmate.bam} \
      -o $BAM_FILE_PATH/${SAMFILE/.sam/_fixmate_sorted.bam}
    
    printf "\n>>> sort done: ${SAMFILE/.sam/_fixmate_sorted.bam} generated"
    rm $BAM_FILE_PATH/${SAMFILE/.sam/_fixmate.bam}
    
    printf "\n>>> initiating index..."

    samtools index \
      -@ 16 \
      $BAM_FILE_PATH/${SAMFILE/.sam/_fixmate_sorted.bam} \
      $BAM_FILE_PATH/${SAMFILE/.sam/_fixmate_sorted.bai}
      
    printf "\n>>> index done: ${SAMFILE/.sam/_fixmate_sorted.bam} indexed\n"
    
  done
```

### 4.1.4 Add read group (@RG) information to BAM files and index them

Fixing for missing Read Group information, a set of reads that were generated from a single run of a sequencing instrument by @RG tag for BAM files. ReadID keyword used in FASTQ files: Illumina HiSeq-HISEQ; MiSeq-M02880 For more Read Groups information, please check https://software.broadinstitute.org/gatk/documentation/article?id=6472

For FASTQ `<file format: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>  <read>:<is filtered>:<control number>:<sample number>>` For more FASTQ file format, please check https://help.basespace.illumina.com/articles/descriptive/fastq-files/

```
RAWFASTQGZ_PATH=$WORKING_PATH/01_SeqData/rawReads.gz
BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/02_BAM
ADD_RG_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/03_addRG_BAM

mkdir -p $ADD_RG_BAM_FILE_PATH

cd $RAWFASTQGZ_PATH

for RAWFASTQGZ_FILE in $(find $RAWFASTQGZ_PATH -name '*R1*.fastq.gz' -exec basename "{}" \;)
  do 
    RGSM_VAR=$(echo $RAWFASTQGZ_FILE | awk -F'_' '{print $1}')
    RGLB_VAR=$(echo $RAWFASTQGZ_FILE | awk -F'_' '{print $1"_LIB001"}')
    RGPU_VAR=$(samtools view -H $RAWFASTQGZ_FILE | grep -Em1 '^@HISEQ|^@M0|^@K' | \
               awk -F':' '{print $3"."$4"."}')$RGSM_VAR
    RGID_VAR=$(samtools view -H $RAWFASTQGZ_FILE | grep -Em1 '^@HISEQ|^@M0|^@K' | awk -F':' '{print $3"."$4"."$2}')
    INPUT_VAR=$BAM_FILE_PATH/$(ls $BAM_FILE_PATH | grep "$RGSM_VAR.*.bam$")
    OUTPUT_VAR=$ADD_RG_BAM_FILE_PATH/$(ls $BAM_FILE_PATH | grep "$RGSM_VAR.*.bam$" | sed 's/.bam/_addRG.bam/g')
    
    gatk AddOrReplaceReadGroups \
      --RGSM $RGSM_VAR     `# {SAMPLE_NAME}` \
      --RGLB $RGLB_VAR     `# {SAMPLE_NAME}_LIB001` \
      --RGPL illumina \
      --RGPU $RGPU_VAR     `# {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_NAME}` \
      --RGID $RGID_VAR     `# {FLOWCELL_BARCODE}.{LANE}.{RUN_NUMBER}` \
      --INPUT $INPUT_VAR \
      --OUTPUT $OUTPUT_VAR
      
    printf '\n>>> initiating index...\n'
    
    samtools index \
    -@ 16 \
    $OUTPUT_VAR \
    ${OUTPUT_VAR/.bam/.bai}
    
    printf ">>> index done: $(find $BAM_FILE_PATH -name "$RGSM_VAR*.bam" | sed 's/.bam/_addRG.bam/g') indexed\n"
    
  done
```

### 4.1.5 Perform local realignment around indels
Local realignment of reads to enhance the alignments in the vicinity of indel polymorphisms
```
ADD_RG_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/03_addRG_BAM
REALIGN_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/04_Realign_BAM
INTERVAL_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/05_Intervals

mkdir -p $REALIGN_BAM_FILE_PATH
mkdir -p $INTERVAL_FILE_PATH

cd $ADD_RG_BAM_FILE_PATH

for ADD_RG_BAM in $(find $ADD_RG_BAM_FILE_PATH -name '*addRG.bam' -exec basename "{}" \;)
  do
    # Identify active intervals for realignments
    gatk3 \
      -T RealignerTargetCreator \
      -R $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      -I $ADD_RG_BAM \
      -o $INTERVAL_FILE_PATH/${ADD_RG_BAM/.bam/.intervals}
    
    # Implement realignments
    gatk3 \
      -T IndelRealigner \
      -R $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      -I $ADD_RG_BAM \
      -targetIntervals $INTERVAL_FILE_PATH/${ADD_RG_BAM/.bam/.intervals} \
      -o $REALIGN_BAM_FILE_PATH/${ADD_RG_BAM/.bam/_realigned.bam}
  done

# Some issues arise from JDK environments problems. Looping will pop-up memory shortage issue and generate 
# damaged or incomplete bam files, instead the previous procedure should implement manually on each sample. 

ERR_LOG_FILE_ARR=($(find . -name '*.log'))

while [ ${#ERR_LOG_FILE_ARR[@]} -gt 0 ]
  do
    for ERR_LOG_FILE in $ERR_LOG_FILE_ARR
      do
    
      # using the error log file to redo previous steps automatically, specifically on the right id of the files.
      # extract id infomation in the error file to redo on the corresponding files
        
        INPUT_FILE=$(grep -o \\-I'\(.*\)'addRG.bam $ERR_LOG_FILE | awk '{print $2}')
        
        gatk3 \
          -T IndelRealigner \
          -R $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
          -I $INPUT_FILE \
          -targetIntervals $INTERVAL_FILE_PATH/${INPUT_FILE/.bam/.intervals} \
          -o $REALIGN_BAM_FILE_PATH/${INPUT_FILE/.bam/_realigned.bam}
      
        rm $ERR_LOG_FILE
      
      done
      
    ERR_LOG_FILE_ARR=($(find . -name '*.log'))
    
  done
```

# 5 Mapping Quality Cotrol
## 5.1 Mark duplicates - using MarkDuplicates from GATK and index them

This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR.

Duplicate reads can also result from a single amplification cluster,incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates. The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).

The tool’s main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024.

This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify and differentiate duplicate types.
```
conda activate mtb_wgs_map_vc_env

WORKING_PATH=[PATH to A WORKING DIREACTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

REF_GENOME_PATH=$REF_PATH/11_Reference_Genome/

REALIGN_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/04_Realign_BAM
DEDUP_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/06_Dedup_BAM

mkdir -p $DEDUP_BAM_FILE_PATH

cd $REALIGN_BAM_FILE_PATH

for REALIGN_BAM in $(find $REALIGN_BAM_FILE_PATH -name "*_realigned.bam" -exec basename "{}" \;)
  do
    BAMFILE_ID=$(echo $REALIGN_BAM | cut -d"_" -f 1)
    
    # picard MarkDuplicates \
    gatk MarkDuplicates \
    --INPUT $REALIGN_BAM \
    --REMOVE_DUPLICATES true \
    --OUTPUT $DEDUP_BAM_FILE_PATH/${REALIGN_BAM/realigned/realigned_dedup} \
    --METRICS_FILE $DEDUP_BAM_FILE_PATH/${BAMFILE_ID}_dedup_metrics.txt
    
    samtools index \
    -@ 16 \
    $DEDUP_BAM_FILE_PATH/${REALIGN_BAM/realigned.bam/realigned_dedup.bam} \
    $DEDUP_BAM_FILE_PATH/${REALIGN_BAM/realigned.bam/realigned_dedup.bai} \
    
    printf "\ndone: ${REALIGN_BAM/realigned.bam/realigned_dedup.bam} indexed\n"
    
  done
```

## 5.2 Base Quality Score Recalibration (BQSR)
First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. The default covariates are read group, reported quality score, machine cycle, and nucleotide context. This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites that are in the known sites VCF. GATK assumes that all reference mismatches we see are therefore errors and indicative of poor base quality. Since there is a large amount of data one can then calculate an empirical probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations.

```
DEDUP_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/06_Dedup_BAM
RECAL_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/07_Recal_BAM
RECAL_TABLE_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/08_Recal_Table
KNOWN_SITE_VCF_PATH=$REF_PATH/12_Calibration_Files

mkdir -p $RECAL_BAM_FILE_PATH $RECAL_TABLE_FILE_PATH

cd $DEDUP_BAM_FILE_PATH

for DEDUP_BAM_FILE in $(find $DEDUP_BAM_FILE_PATH -name "*.bam" -exec basename "{}" \;)
  do
    
    gatk BaseRecalibrator \
      --input $DEDUP_BAM_FILE \
      --known-sites $(find $KNOWN_SITE_VCF_PATH -name "*_mdf.vcf") \
      --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      --output $RECAL_TABLE_FILE_PATH/${DEDUP_BAM_FILE/.bam/_recal.table}
    
    gatk ApplyBQSR \
      --input $DEDUP_BAM_FILE \
      --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      --bqsr-recal-file $RECAL_TABLE_FILE_PATH/${DEDUP_BAM_FILE/.bam/_recal.table} \
      --output $RECAL_BAM_FILE_PATH/${DEDUP_BAM_FILE/.bam/_recal.bam}
    
    gatk BaseRecalibrator \
      --input $RECAL_BAM_FILE_PATH/${DEDUP_BAM_FILE/.bam/_recal.bam} \
      --known-sites $(find $KNOWN_SITE_VCF_PATH -name "*_mdf.vcf") \
      --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      --output $RECAL_TABLE_FILE_PATH/${DEDUP_BAM_FILE/.bam/_post_recal.table}
    
    gatk ApplyBQSR \
      --input $RECAL_BAM_FILE_PATH/${DEDUP_BAM_FILE/.bam/_recal.bam} \
      --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      --bqsr-recal-file $RECAL_TABLE_FILE_PATH/${DEDUP_BAM_FILE/.bam/_post_recal.table} \
      --output $RECAL_BAM_FILE_PATH/${DEDUP_BAM_FILE/.bam/_post_recal.bam}
    
  done

```

## 5.3 Alignment quality profiling by Qualimap
```
conda activate mtb_wgs_qc_env

RECAL_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/07_Recal_BAM
BAM_QC_OUTPUT_PATH=$WORKING_PATH/02_Quality_Control/bam_qc_report
BAM_MULTI_QC_OUTPUT_PATH=$WORKING_PATH/02_Quality_Control/bam_multi_qc_report

mkdir -p $BAM_QC_OUTPUT_PATH $BAM_MULTI_QC_OUTPUT_PATH

# Individual BAM QC reports
cd $RECAL_BAM_FILE_PATH

for RECAL_BAM in $(find $RECAL_BAM_FILE_PATH -name "*_dedup_recal.bam" -exec basename "{}" \;)
  do
    qualimap bamqc \
      -bam $RECAL_BAM \
      -outdir $BAM_QC_OUTPUT_PATH/$RECAL_BAM \
      --java-mem-size=64G
  done

# Multiple BAM QC report
cd $RECAL_BAM_FILE_PATH

# Sample file path column
paste -d'\t' \
<(find $RECAL_BAM_FILE_PATH -name "*_dedup_recal.bam" -exec basename "{}" \; | cut -d'_' -f 1) `# Sample name` \
<(find $BAM_QC_OUTPUT_PATH -name "*_dedup_recal.bam") `# Sample file path` \
> $BAM_MULTI_QC_OUTPUT_PATH/multi_bam_qc_report_input.txt

qualimap multi-bamqc \
  --data $BAM_MULTI_QC_OUTPUT_PATH/multi_bam_qc_report_input.txt \
  --feature-file $(find $REF_GENOME_PATH -maxdepth 5 -name "*.gff") \
  -outdir $BAM_MULTI_QC_OUTPUT_PATH \
  --java-mem-size=64G
```

# 6. Variant Calling
Available tools: samtools(bcftools), GATK-HaplotypeCaller, CRISP, VarScan2

## 6.1 Variant calling by samtools (bcftools)
Joint alignments for variant calling
```
conda activate mtb_wgs_map_vc_env

WORKING_PATH=[PATH to A WORKING DIRECTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

REF_GENOME_PATH=$REF_PATH/11_Reference_Genome/

RECAL_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/07_Recal_BAM
BCFTOOLS_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/01_bcftools
VARIANT_QC_REPORT_PATH=$WORKING_PATH/02_Quality_Control/variant_qc_report

mkdir -p $BCFTOOLS_VCF_FILE_PATH $VARIANT_QC_REPORT_PATH

cd $RECAL_BAM_FILE_PATH

bcftools mpileup \
  $(find $RECAL_BAM_FILE_PATH -name "*_dedup_recal.bam" -exec basename "{}" \;) \
  --output-type u \
  --fasta-ref $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --threads 16 |
bcftools call \
  --ploidy 1 \
  --multiallelic-caller  `# alternative model for multiallelic and rare-variant calling` \
  --variants-only        `# output variant sites only` \
  --output-type z        `# output compressed VCF` \
  --output $BCFTOOLS_VCF_FILE_PATH/whole_var.vcf.gz \
  --threads 16


tabix \
  --force \
  --preset vcf \
  $BCFTOOLS_VCF_FILE_PATH/whole_var.vcf.gz

bcftools stats \
  $BCFTOOLS_VCF_FILE_PATH/whole_var.vcf.gz \
  --fasta-ref $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --samples - \
  --threads 16 \
  > $VARIANT_QC_REPORT_PATH/samtools_vcf.stats

plot-vcfstats \
  $VARIANT_QC_REPORT_PATH/samtools_vcf.stats \
  --prefix $VARIANT_QC_REPORT_PATH/samtools_vcf_stats_plot/
```

## 6.2 Variant calling by GATK-HaplotypeCaller
Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect.

### 6.2.1 Creating the fasta sequence dictionary file
```
REF_GENOME_PATH=$REF_PATH/11_Reference_Genome/

cd $REF_GENOME_PATH

gatk CreateSequenceDictionary \
  --REFERENCE $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --OUTPUT $(find $REF_GENOME_PATH -name "*.genomic.fna" | sed 's/.fna/.dict/')
```

### 6.2.2 GATK gVCF mode for variant calling individually
```
RECAL_BAM_FILE_PATH=$WORKING_PATH/03_Seq_Alignment_Map/01_BWA_MEM/07_Recal_BAM
GATK_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/02_GATK/01_haplotypecaller_gvcf

mkdir -p $GATK_VCF_FILE_PATH

cd $RECAL_BAM_FILE_PATH

for RECAL_BAM_FILE in $(find $RECAL_BAM_FILE_PATH -name "*_dedup_recal.bam" -exec basename "{}" \;)
  do
    OUTPUT=${RECAL_BAM_FILE/fixmate_sorted_addRG_realigned_dedup_recal.bam/raw.snps.indels.g.vcf}
    
    gatk HaplotypeCaller \
      --input $RECAL_BAM_FILE \
      --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      --emit-ref-confidence GVCF \
      --read-filter NotOpticalDuplicateReadFilter \
      --sample-ploidy 1 \
      --output $GATK_VCF_FILE_PATH/$OUTPUT
  done
```

### 6.2.3 Consolidate GVCFs
```
GVCF_GENOMICS_DB_WORK_PATH=$WORKING_PATH/04_Variant_Calling/02_GATK/02_gvcf_geno_db

cd $GATK_VCF_FILE_PATH

[ $(find .. $GATK_VCF_FILE_PATH -d -name '02_gvcf_geno_db') ] && rm -r "../02_gvcf_geno_db"

gatk GenomicsDBImport \
  $(find $GATK_VCF_FILE_PATH -name "*.g.vcf" -exec basename "{}" \; | awk '{print "--variant"" "$1" "}') \
  --genomicsdb-workspace-path $GVCF_GENOMICS_DB_WORK_PATH \
  --intervals "$(find $REF_GENOME_PATH -name "*.genomic.fna" -exec basename "{}" \; | cut -f1,2 -d'_')" \
  --reader-threads 16 \
  --overwrite-existing-genomicsdb-workspace
```

### 6.2.4 Joint-call cohort
```
GATK_JOINT_CALL_COHORT_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/02_GATK/03_joint_call_cohort_vcf

mkdir -p $GATK_JOINT_CALL_COHORT_VCF_FILE_PATH

cd $GVCF_GENOMICS_DB_WORK_PATH

gatk GenotypeGVCFs \
  --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --variant gendb://$GVCF_GENOMICS_DB_WORK_PATH \
  --use-new-qual-calculator \
  --output $GATK_JOINT_CALL_COHORT_VCF_FILE_PATH/'joint_call_cohort.vcf.gz'
```

# 7. Variants Filtration
## 7.1 (Currently unable to implement …) Variant Quality Score Recalibration (VQSR)
```
Currently no available variant database to train model. (2019.10.20)
```

## 7.2 Extract common variants between VCF files generated from BCFtools and GATK
```
conda activate mtb_wgs_map_vc_env

WORKING_PATH=[PATH to A WORKING DIRECTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

REF_GENOME_PATH=$REF_PATH/11_Reference_Genome/

BCFTOOLS_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/01_bcftools/
GATK_JOINT_CALL_COHORT_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/02_GATK/03_joint_call_cohort_vcf
COMMON_VARIANT_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/00_common_variant_vcf

mkdir -p $COMMON_VARIANT_VCF_FILE_PATH

bcftools isec \
 $(find $BCFTOOLS_VCF_FILE_PATH -name "*vcf.gz") \
 $(find $GATK_JOINT_CALL_COHORT_VCF_FILE_PATH -name "*vcf.gz") \
 --threads 16 \
 --prefix $COMMON_VARIANT_VCF_FILE_PATH \
 --output-type z  # compressed VCF
```

## 7.3 Variant calling quality profiling for raw variants by DISCVRseq
```
alias DISCVRseq='java -jar $TOOLS_PATH/DISCVR-seq_Toolkit/DISCVRSeq-1.07.jar'

BCFTOOLS_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/01_bcftools/
GATK_JOINT_CALL_COHORT_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/02_GATK/03_joint_call_cohort_vcf
COMMON_VARIANT_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/00_common_variant_vcf
VARIANT_QC_REPORT_PATH=$WORKING_PATH/02_Quality_Control/variant_qc_report

mkdir -p $VARIANT_QC_REPORT_PATH

cd $VARIANT_QC_REPORT_PATH

DISCVRseq VariantQC \
  --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --variant $(find $BCFTOOLS_VCF_FILE_PATH -name "*.vcf.gz") \
  --output $VARIANT_QC_REPORT_PATH/00_samtools_vcf_qc.html

DISCVRseq VariantQC \
  --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --variant $(find $GATK_JOINT_CALL_COHORT_VCF_FILE_PATH -name "*.vcf.gz") \
  --output $VARIANT_QC_REPORT_PATH/01_GATK_gvcf_qc.html

DISCVRseq VariantQC \
  --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --variant $(find $COMMON_VARIANT_VCF_FILE_PATH -name "0002.vcf.gz") \
  --output $VARIANT_QC_REPORT_PATH/02_samtools_GATK_common_vcf_qc.html

DISCVRseq VariantQC \
  --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --variant $(find $COMMON_VARIANT_VCF_FILE_PATH -name "0003.vcf.gz") \
  --output $VARIANT_QC_REPORT_PATH/03_GATK_samtools_common_vcf_qc.html
```

## 7.4 Raw variant set evaluation
### 7.4.1 Raw variants by PICARD and GATK
```
KNOWN_SITE_VCF_PATH=$REF_PATH/12_Calibration_Files

COMMON_VARIANT_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/00_common_variant_vcf
RAW_VARIANT_EVAL_PATH=$WORKING_PATH/02_Quality_Control/raw_variant_evaluation

mkdir -p $RAW_VARIANT_EVAL_PATH

cd $RAW_VARIANT_EVAL_PATH

picard CollectVariantCallingMetrics \
  INPUT=$COMMON_VARIANT_VCF_FILE_PATH/0003.vcf.gz \
  OUTPUT=$RAW_VARIANT_EVAL_PATH/GATK_bcftools_Metrics_byPICARD.raw \
  DBSNP=$(find $KNOWN_SITE_VCF_PATH -name "MTB_Base_Calibration_List_mdf.vcf")

gatk VariantEval \
  --reference $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
  --eval $COMMON_VARIANT_VCF_FILE_PATH/0003.vcf.gz \
  --comp $(find $KNOWN_SITE_VCF_PATH -name "MTB_Base_Calibration_List_mdf.vcf") `# inappropriate comparison` \
  --do-not-use-all-standard-modules \
  --eval-module CompOverlap \
  --eval-module IndelSummary \
  --eval-module TiTvVariantEvaluator \
  --eval-module CountVariants -EV MultiallelicSummary \
  --output $RAW_VARIANT_EVAL_PATH/GATK_bcftools_Metrics_byGATK.raw_variant_calling_metrics
```

### 7.4.2 Visualize raw variants variant calling evaluation of PICARD by R
#### 7.4.2.1 `[R]`data import
```
.libPaths(c("~/3_Software/3_02_Programming/R_Packages/bioc_lib/",
            Sys.getenv("R_LIBS_USER")))

library(data.table)
library(magrittr)
library(stringr)
library(ggplot2)
library(VariantAnnotation)

WORKING_PATH          <- [PATH to A WORKING DIRECTORY]
REF_PATH              <- paste0(str_extract(WORKING_PATH,pattern = ".*.tuberculosis/"),"00_Reference_File/11_Reference_Genome")

RAW_VARIANT_EVAL_PATH <- paste0(WORKING_PATH,"02_Quality_Control/raw_variant_evaluation")
COMMON_VCF_FILE_PATH  <- paste0(WORKING_PATH,"04_Variant_Calling/00_common_variant_vcf")

common_summary_mtc <- fread(list.files(RAW_VARIANT_EVAL_PATH,pattern = "summary_metrics",full.names = TRUE))
common_detail_mtc  <- fread(list.files(RAW_VARIANT_EVAL_PATH,pattern = "detail_metrics",full.names = TRUE))
```
#### 7.4.2.2 `[R]`data manipulation
```
common_VCF_annotation_value <- 
  readVcf(file = list.files(path = COMMON_VCF_FILE_PATH,pattern = "0003.vcf.gz$",full.names = TRUE),
          genome = list.files(REF_PATH,pattern = "*.gff$",recursive = TRUE,full.names = TRUE)) %>% 
  { as.data.table(.@info) } %>% 
  .[,.SD,.SDcols = c("AF","MLEAF","QD","FS","SOR","MQ","MQRankSum","ReadPosRankSum")] %>% 
  .[,':=' (AF = unlist(AF),MLEAF = unlist(MLEAF))] %>% 
  melt(.,measure.vars = c("AF","MLEAF","QD","FS","SOR","MQ","MQRankSum","ReadPosRankSum")) %>% 
  .[,variable := factor(variable,labels = c("AlleleFrequency","MLEAlleleFrequency",
                                            "QualByDepth","FisherStrand","StraindOddRatio","RMSMappingQuality",
                                            "MappingQualityRankSumTest","ReadPosRankSumTest"))]
```
#### 7.4.2.3 `[R]`data visualization
```
ggplot(data = common_VCF_annotation_value,aes(x = value,y = stat(scaled))) +
  geom_density(na.rm = TRUE) +
  facet_wrap( ~ common_VCF_annotation_value$variable,
              scales = "free",nrow = 3,ncol = 3)
```

### 7.5 Hard filtering for SNP and INDEL from common variants VCF file
```
COMMON_VARIANT_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/00_common_variant_vcf
FILTERED_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/90_filtered_vcf

mkdir -p $FILTERED_VCF_FILE_PATH

cd $FILTERED_VCF_FILE_PATH

# SNP filtering
gatk SelectVariants \
  --variant $(find $COMMON_VARIANT_VCF_FILE_PATH -name "0003.vcf.gz") \
  --select-type-to-include SNP \
  --output $FILTERED_VCF_FILE_PATH/common_SNP.vcf.gz

gatk VariantFiltration \
  --variant $(find $FILTERED_VCF_FILE_PATH -name "common_SNP.vcf.gz") \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --output $FILTERED_VCF_FILE_PATH/common_SNP_filter_status_annotated.vcf.gz

gatk SelectVariants \
  --variant $(find $FILTERED_VCF_FILE_PATH -name "*SNP_filter_status_annotated.vcf.gz") \
  --selectExpressions 'vc.isNotFiltered()' \
  --output $FILTERED_VCF_FILE_PATH/common_SNP_hard_filtered.vcf.gz

# INDEL filtering
gatk SelectVariants \
  --variant $(find $COMMON_VARIANT_VCF_FILE_PATH -name "0003.vcf.gz") \
  --select-type-to-include INDEL \
  --output $FILTERED_VCF_FILE_PATH/common_INDEL.vcf.gz

gatk VariantFiltration \
  --variant $(find $FILTERED_VCF_FILE_PATH -name "common_INDEL.vcf.gz") \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  --output $FILTERED_VCF_FILE_PATH/common_INDEL_filter_status_annotated.vcf.gz

gatk SelectVariants \
  --variant $(find $FILTERED_VCF_FILE_PATH -name "*INDEL_filter_status_annotated.vcf.gz") \
  --selectExpressions 'vc.isNotFiltered()' \
  --output $FILTERED_VCF_FILE_PATH/common_INDEL_hard_filtered.vcf.gz
```

### 7.6 Venn diagrams for raw and filtered SNPs and INDELs called by GATK and samtools
#### 7.6.1 Generate position list from each VCF files for venn diagram
```
BCFTOOLS_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/01_bcftools/
GATK_JOINT_CALL_COHORT_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/02_GATK/03_joint_call_cohort_vcf
FILTERED_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/90_filtered_vcf
VCF_POSLIST_FILE_PATH=$WORKING_PATH/04_Variant_Calling/99_vcf_poslist

mkdir -p $VCF_POSLIST_FILE_PATH

cd $VCF_POSLIST_FILE_PATH

# samtools SNP position list extraction
bcftools view \
  --type snps \
  $(find $BCFTOOLS_VCF_FILE_PATH -name "whole_var.vcf.gz") \
  | awk '/^[^#]/ {print $2}' \
  > $VCF_POSLIST_FILE_PATH/1_samtools_raw_SNP_poslist.txt

# GATK SNP position list extraction
bcftools view \
  --type snps \
  $(find $GATK_JOINT_CALL_COHORT_VCF_FILE_PATH -name "joint_call_cohort.vcf.gz") \
  | awk '/^[^#]/ {print $2}' \
  > $VCF_POSLIST_FILE_PATH/2_GATK_raw_SNP_poslist.txt

# GATK-samtools common and filtered SNP position list extraction
bcftools view \
  --type snps \
  $(find $FILTERED_VCF_FILE_PATH -name "*SNP_hard_filtered.vcf.gz") \
  | awk '/^[^#]/ {print $2}' \
  > $VCF_POSLIST_FILE_PATH/3_GATK_samtools_hard_filtered_SNP_poslist.txt

# samtools INDEL list extraction
bcftools view \
  --type indels \
  $(find $BCFTOOLS_VCF_FILE_PATH -name "whole_var.vcf.gz") \
  | awk '/^[^#]/ {print $2}' \
  > $VCF_POSLIST_FILE_PATH/4_samtools_raw_INDEL_poslist.txt

# GATK SNP list extraction
bcftools view \
  --type indels \
  $(find $GATK_JOINT_CALL_COHORT_VCF_FILE_PATH -name "joint_call_cohort.vcf.gz") \
  | awk '/^[^#]/ {print $2}' \
  > $VCF_POSLIST_FILE_PATH/5_GATK_raw_INDEL_poslist.txt

# GATK-samtools common and filtered INDEL position list extraction
bcftools view \
  --type indels \
  $(find $FILTERED_VCF_FILE_PATH -name "*INDEL_hard_filtered.vcf.gz") \
  | awk '/^[^#]/ {print $2}' \
  > $VCF_POSLIST_FILE_PATH/6_GATK_samtools_hard_filtered_INDEL_poslist.txt
```

### 7.6.2 `[R]`Plot Venn diagrams for SNP sets
```
library(data.table)
library(magrittr)
library(stringr)
library(ggplot2)
library(VennDiagram)
library(png)
library(cowplot)

WORKING_PATH          <- [PATH to A WORKING DIRECTORY]
VCF_POSLIST_FILE_PATH <- paste0(WORKING_PATH,"04_Variant_Calling/99_vcf_poslist")

Category <- c("SAMtools","GATK","Filtered common set")
SNPset <- list.files(VCF_POSLIST_FILE_PATH,pattern = "*SNP_poslist.txt",full.names = TRUE) %>% 
          lapply(.,function(x) fread(x) %>% unlist(.,use.names = FALSE))

# Venn diagrams
venn.diagram(x = list(SNPset[[1]],SNPset[[2]],SNPset[[3]]),
             category.names = Category,
             
             # Circles
             col=c("#440154ff","#21908dff","#fde725ff"),
             fill = c(alpha("#440154ff",0.4), alpha("#21908dff",0.4), alpha("#fde725ff",0.4)),
             euler.d = FALSE,
             scaled = FALSE,
             
             # Numbers
             cex = rel(1.3),
             fontface = "bold",
             fontfamily = "Arial",
             
             # Set names
             cat.cex = rel(1.3),
             cat.fontface = "bold",
             cat.fonfamily = "Arial",
             cat.just = list(c(0.5,1),c(0.5,1.5),c(0.5,0.5)),
             
             # Output features
             filename = paste(VCF_POSLIST_FILE_PATH,"SNPset_VennDiagram.png",sep = "/"),
             imagetype = "png")

list.files(VCF_POSLIST_FILE_PATH,"*.log",full.names = TRUE) %>% file.remove(.)

venn.diagram.png <- readPNG(list.files(VCF_POSLIST_FILE_PATH,"SNPset_VennDiagram.png$",full.names = TRUE))

# Barpolot
count_data <- data.table(category = factor(Category,levels = Category),
                         count = sapply(SNPset,function(x) length(x)))

barplot <- 
  ggplot(data = count_data) +
    geom_bar(aes(x = category,y = count,fill = category),stat = "identity",alpha = 0.5) +
    geom_text(aes(x = category,y = count,label = count,vjust = -5),
              color = "black",fontface = "bold",size = rel(3),position = position_stack(vjust = 0.85)) +
    scale_fill_manual(values = c("#440154ff", '#21908dff', '#fde725ff')) +
    scale_x_discrete(name = "Variant call sets") +
    scale_y_continuous(name = "Variant counts") +
    coord_flip() +
    theme(panel.background = element_blank(),
          axis.line = element_line(),
          axis.title.x = element_text(face = "bold",size = rel(.9),margin = margin(t = 10,r = 0,b = 0,l = 0)),
          axis.title.y = element_text(face = "bold",size = rel(.9),margin = margin(t = 0,r = 10,b = 0,l = 0)),
          axis.text = element_text(face = "bold",size = rel(.8)),
          legend.position = "none")

plot_grid(rasterGrob(venn.diagram.png),barplot,rel_heights = c(3,1),ncol = 2,
          labels = c("Venn diagram of SNP counts","SNP counts"),hjust = -0.02,vjust = 1.5)
```

### 7.6.3 `[R]`Plot Venn diagrams for INDEL sets
```
library(data.table)
library(magrittr)
library(ggplot2)
library(VennDiagram)
library(png)
library(cowplot)

WORKING_PATH          <- [PATH to A WORKING DIRECTORY]
VCF_POSLIST_FILE_PATH <- paste0(WORKING_PATH,"04_Variant_Calling/99_vcf_poslist")

Category <- c("SAMtools","GATK","Filtered common set")
INDELset <- list.files(VCF_POSLIST_FILE_PATH,pattern = "*INDEL_poslist.txt",full.names = TRUE) %>% 
            lapply(.,function(x) fread(x) %>% unlist(.,use.names = FALSE))

# Venn diagrams
venn.diagram(x = list(INDELset[[1]],INDELset[[2]],INDELset[[3]]),
             category.names = Category,
             
             # Circles
             col=c("#B3E2CD","#FDCDAC","#CBD5E8"),
             fill = c(alpha("#B3E2CD",0.8), alpha("#FDCDAC",0.8), alpha("#CBD5E8",0.8)),
             euler.d = FALSE,
             scaled = FALSE,
             
             # Numbers
             cex = rel(1.3),
             fontface = "bold",
             fontfamily = "Arial",
             
             # Set names
             cat.cex = rel(1.3),
             cat.fontface = "bold",
             cat.fonfamily = "Arial",
             cat.just = list(c(0.5,1),c(0.5,1.5),c(0.5,0.5)),
             
             # Output features
             filename = paste(VCF_POSLIST_FILE_PATH,"INDELset_VennDiagram.png",sep = "/"),
             imagetype = "png")

list.files(VCF_POSLIST_FILE_PATH,"*.log",full.names = TRUE) %>% file.remove(.)

venn.diagram.png <- readPNG(list.files(VCF_POSLIST_FILE_PATH,"INDELset_VennDiagram.png$",full.names = TRUE))

# Barpolot
count_data <- data.table(category = factor(Category,levels = Category),
                         count = sapply(INDELset,function(x) length(x)))

barplot <- 
  ggplot(data = count_data) +
    geom_bar(aes(x = category,y = count,fill = category),stat = "identity",alpha = 0.9) +
    geom_text(aes(x = category,y = count,label = count,vjust = -5),
              color = "black",fontface = "bold",size = rel(3),position = position_stack(vjust = 0.85)) +
    scale_fill_manual(values = c("#B3E2CD","#FDCDAC","#CBD5E8")) +
    scale_x_discrete(name = "Variant call sets") +
    scale_y_continuous(name = "Variant counts") +
    coord_flip() +
    theme(panel.background = element_blank(),
          axis.line = element_line(),
          axis.title.x = element_text(face = "bold",size = rel(.9),margin = margin(t = 10,r = 0,b = 0,l = 0)),
          axis.title.y = element_text(face = "bold",size = rel(.9),margin = margin(t = 0,r = 10,b = 0,l = 0)),
          axis.text = element_text(face = "bold",size = rel(.8)),
          legend.position = "none")

plot_grid(rasterGrob(venn.diagram.png),barplot,rel_heights = c(3,1),ncol = 2,
          labels = c("Venn diagram of INDEL counts","INDEL counts"),hjust = -0.02,vjust = 1.5)
```

## 7.7 Filtered variant set evaluation
### 7.7.1 Filtered variants by PICARD and GATK
```
FILTERED_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/90_filtered_vcf
FILTERED_VARIANT_EVAL_PATH=$WORKING_PATH/02_Quality_Control/filtered_variant_evaluation

REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/
REF_GENOME_PATH=$REF_PATH/11_Reference_Genome/
KNOWN_SITE_VCF_PATH=$REF_PATH/12_Calibration_Files


mkdir -p $FILTERED_VARIANT_EVAL_PATH

cd $FILTERED_VARIANT_EVAL_PATH

picard CollectVariantCallingMetrics \
  INPUT=$(find $FILTERED_VCF_FILE_PATH -name "*SNP_hard_filtered.vcf.gz") \
  OUTPUT=$FILTERED_VARIANT_EVAL_PATH/GATK_samtools_Metrics_byPICARD.filtered.SNP \
  DBSNP=$(find $KNOWN_SITE_VCF_PATH -name "MTB_Base_Calibration_List_mdf.vcf")

gatk VariantEval \
  --reference $(find $REF_GENOME_PATH -name "*.genomic.fna") \
  --eval $(find $FILTERED_VCF_FILE_PATH -name "*SNP_hard_filtered.vcf.gz") \
  --comp $(find $KNOWN_SITE_VCF_PATH -name "MTB_Base_Calibration_List_mdf.vcf") `# inappropriate comparison` \
  --do-not-use-all-standard-modules \
  --eval-module CompOverlap \
  --eval-module IndelSummary \
  --eval-module TiTvVariantEvaluator \
  --eval-module CountVariants -EV MultiallelicSummary \
  --output $FILTERED_VARIANT_EVAL_PATH/GATK_samtools_Metrics_byGATK.filtered.SNP.variant_calling_metrics
```

### 7.7.2 `[R]`Visualize filtered variants variant calling evaluation of PICARD by R
#### 7.7.2.1 `[R]`data import
```
.libPaths(c("~/3_Software/3_02_Programming/R_Packages/bioc_lib/",
            Sys.getenv("R_LIBS_USER")))

library(data.table)
library(magrittr)
library(ggplot2)
library(VariantAnnotation)

WORKING_PATH               <- [PATH to A WORKING DIRECTORY]
FILTERED_VARIANT_EVAL_PATH <- paste0(WORKING_PATH,"02_Quality_Control/filtered_variant_evaluation")
FILTERED_VCF_FILE_PATH     <- paste0(WORKING_PATH,"04_Variant_Calling/90_filtered_vcf")


summary_mtc <- fread(list.files(FILTERED_VARIANT_EVAL_PATH,pattern = "summary_metrics",full.names = TRUE))
detail_mtc  <- fread(list.files(FILTERED_VARIANT_EVAL_PATH,pattern = "detail_metrics",full.names = TRUE))
```

#### 7.7.2.2 `[R]`data manipulation for SNP and INDEL call set
```
common_SNP_VCF_annotation_value <- 
  readVcf(file = list.files(path = FILTERED_VCF_FILE_PATH,
                            pattern = "SNP_hard_filtered.vcf.gz$",full.names = TRUE),
          genome = paste0(WORKING_PATH,
                          "09_Reference_Genome/refseq/bacteria/GCF_000195955.2/",
                          "NC_000962.3_Mtb_H37Rv.genomic.gff")) %>% 
  { as.data.table(.@info) } %>% 
  .[,.SD,.SDcols = c("QD","FS","SOR","MQ","MQRankSum","ReadPosRankSum")] %>% 
  # .[,':=' (AF = unlist(AF),MLEAF = unlist(MLEAF))] %>% 
  melt(.,measure.vars = c("QD","FS","SOR","MQ","MQRankSum","ReadPosRankSum")) %>% 
  .[,variable := factor(variable,labels = c("QualByDepth","FisherStrand","StraindOddRatio","RMSMappingQuality",
                                            "MappingQualityRankSumTest","ReadPosRankSumTest"))]

common_INDEL_VCF_annotation_value <- 
  readVcf(file = list.files(path = FILTERED_VCF_FILE_PATH,
                            pattern = "INDEL_hard_filtered.vcf.gz$",full.names = TRUE),
          genome = paste0(WORKING_PATH,
                          "09_Reference_Genome/refseq/bacteria/GCF_000195955.2/",
                          "NC_000962.3_Mtb_H37Rv.genomic.gff")) %>% 
  { as.data.table(.@info) } %>% 
  .[,.SD,.SDcols = c("QD","FS","ReadPosRankSum")] %>% 
  melt(.,measure.vars = c("QD","FS","ReadPosRankSum")) %>% 
  .[,variable := factor(variable,labels = c("QualByDepth","FisherStrand","ReadPosRankSumTest"))]
```

#### 7.7.2.3 `[R]`data visualization for SNP call set
```
ggplot(data = common_SNP_VCF_annotation_value,aes(x = value,y = stat(scaled))) +
  geom_density(na.rm = TRUE) +
  scale_y_continuous(name = "Scaled density") +
  facet_wrap( ~ common_SNP_VCF_annotation_value$variable,
              scales = "free",ncol = 3)

# ggsave(filename = "~/Desktop/VCF_annotation_value_profile_simple.jpeg",
#        device = "jpeg",width = 192,height = 168,units = "mm",dpi = 600)
```

#### 7.7.2.4 `[R]`data visualization for INDEL call set
```
ggplot(data = common_INDEL_VCF_annotation_value,aes(x = value,y = stat(scaled))) +
  geom_density(na.rm = TRUE) +
  scale_y_continuous(name = "Scaled density") +
  facet_wrap( ~ common_INDEL_VCF_annotation_value$variable,
              scales = "free",ncol = 3,drop = FALSE)
```

# 8. Variants Annotation by SnpEff
```
conda activate mtb_wgs_anno_env

WORKING_PATH=[PATH to A WORKING DIRECTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

GFF_FILE_PATH=$REF_PATH/11_Reference_Genome/
FILTERED_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/90_filtered_vcf

cd $FILTERED_VCF_FILE_PATH

# Prepare the annotation database
snpeff databases | awk -F'\t' '{print $1"\t"$2"\t"$3}' | grep Mycobacterium_tuberculosis_h37rv | less
VARIANTS_DATABASE='Mycobacterium_tuberculosis_h37rv'
snpeff download -v $(echo $VARIANTS_DATABASE)

# Convert the column name and its values for the issue of index mapping
bcftools view ./common_SNP_hard_filtered.vcf.gz | \
 sed 's/NC_000962.3/Chromosome/g' \
 > ./converted_common_SNP_hard_filtered.vcf

# Initiate the annotation procedure
snpeff ann \
  $(echo $VARIANTS_DATABASE)                          `# genome version` \
  ./converted_common_SNP_hard_filtered.vcf            `# input file`\
  -noLof                                              `# Do not add LOF and NMD annotations` \
  -upDownStreamLen 0                                  `# Set upstream downstream interval length (in bases)` \
  -no-downstream                                      `# Do not show DOWNSTREAM changes` \
  -noLog                                              `# Do not report usage statistics to server` \
 > ./converted_common_SNP_hard_filtered_annotated.vcf `# output to a new vcf`

# Revert the column name and its values
cat ./converted_common_SNP_hard_filtered_annotated.vcf | \
  sed -e 's/Chromosome/NC_000962.3/g' \
  > common_SNP_hard_filtered_annotated.vcf

# Compress the vcf
bcftools view \
  --compression-level 4 \
  --output-type z \
  --output-file common_SNP_hard_filtered_annotated.vcf.gz \
  common_SNP_hard_filtered_annotated.vcf

# Index the vcf  
tabix \
  --force \
  --preset vcf \
  common_SNP_hard_filtered_annotated.vcf.gz

```

# 9. Filter Variants by Gene or Function
## 9.1 Filter variants within PE/PPE region
```
conda activate mtb_wgs_anno_env

WORKING_PATH=[PATH to A WORKING DIRECTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

FILTERED_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/90_filtered_vcf
FILTERED_BYGENE_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/91_filtered_byGene_vcf

mkdir -p $FILTERED_BYGENE_VCF_FILE_PATH

cd $FILTERED_VCF_FILE_PATH

# Counts of variant on PE/PPE regions
cat ./common_SNP_hard_filtered_annotated.vcf | \
  grep -o '.*.ANN=.*.PE.*.GT\:AD' | \
  wc -l

# Remove variants on PE/PPE regions
cat ./common_SNP_hard_filtered_annotated.vcf | \
  snpsift filter -n "( ANN[*].GENE =~ 'PE' )"  | \
  > $FILTERED_BYGENE_VCF_FILE_PATH/01_rm_PE-PPE_common_SNP_annotated.vcf
  
# Compress the vcf
bcftools view \
  --compression-level 4 \
  --output-type z \
  --output-file $FILTERED_BYGENE_VCF_FILE_PATH/01_rm_PE-PPE_common_SNP_annotated.vcf.gz \
  $FILTERED_BYGENE_VCF_FILE_PATH/01_rm_PE-PPE_common_SNP_annotated.vcf

# Index the vcf  
tabix \
  --force \
  --preset vcf \
  $FILTERED_BYGENE_VCF_FILE_PATH/01_rm_PE-PPE_common_SNP_annotated.vcf.gz
```

## 9.2 Filter variants pertaining to drug resistance
```
conda activate mtb_wgs_anno_env

WORKING_PATH=[PATH to A WORKING DIRECTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

VARIANT_LIST_FILE=$(find $REF_PATH/12_Calibration_Files -name "**Base_Calibration_List_rmPhylo.vcf.gz")

FILTERED_BYGENE_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/91_filtered_byGene_vcf

# Remove variants pertaining to drug resistance
bcftools isec \
 $(find $FILTERED_BYGENE_VCF_FILE_PATH -name "01_*vcf.gz") \
 $VARIANT_LIST_FILE \
 --complement \
 --threads 16 \
 --prefix $FILTERED_BYGENE_VCF_FILE_PATH \
 --output-type z  # compressed VCF

# Index the vcf file
tabix \
  --force \
  --preset vcf \
  0000.vcf.gz

# Rename the vcf file
cd $FILTERED_BYGENE_VCF_FILE_PATH
mv 0000.vcf.gz 02_rm_DRG_common_SNP_annotated.vcf.gz
mv 0000.vcf.gz.tbi 02_rm_DRG_common_SNP_annotated.vcf.gz.tbi
```

# 10. Generate Concensus Sequence
## 10.1 FASTA files
For each samples based on VCF file against reference genome sequence
```
conda activate mtb_wgs_anno_env

WORKING_PATH=[PATH to A WOKRING DIRECTORY]
REF_PATH=$(echo $WORKING_PATH | awk -F'/tuberculosis' '{print $1"/tuberculosis"}')/00_Reference_File/

REF_GENOME_PATH=$REF_PATH/11_Reference_Genome/

# [!RUN EITHER] IF APPLY Unfiltered-by-Gene VCF FILE
FILTERED_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/90_filtered_vcf/
TARGET_VCF_FILE=$(find $FILTERED_VCF_FILE_PATH -name "*common_SNP_annotated.vcf")
CONSENSUS_FASTA_FILE=$WORKING_PATH/05_Concensus_Fasta

# [!RUN EITHER] IF APPLY Filtered-by-Gene VCF FILE
FILTERED_BYGENE_VCF_FILE_PATH=$WORKING_PATH/04_Variant_Calling/91_filtered_byGene_vcf/
TARGET_VCF_FILE=$(find $FILTERED_BYGENE_VCF_FILE_PATH -name "02*common_SNP_annotated.vcf.gz")
CONSENSUS_FASTA_FILE=$WORKING_PATH/05_Concensus_Fasta

mkdir -p $CONSENSUS_FASTA_FILE

cd $FILTERED_BYGENE_VCF_FILE_PATH

# Generate consensus sequence for each samples
for SAMPLE_ID in $(bcftools query -l $TARGET_VCF_FILE)
  do 
    bcftools consensus \
      --fasta-ref $(find $REF_GENOME_PATH -maxdepth 5 -name "*Mtb_H37Rv.genomic.fna") \
      --haplotype A \
      --sample $SAMPLE_ID \
      --prefix "SampleID:$SAMPLE_ID consensus sequence against " \
      --output $CONSENSUS_FASTA_FILE/$SAMPLE_ID.fasta \
      $TARGET_VCF_FILE
  done

# Concatenate sequence data into one fasta and calculate distances
cd $CONSENSUS_FASTA_FILE

cat *.fasta > merged_seq.fasta
snp-dists -b merged_seq.fasta > merged_seq_distances.txt

```

## 10.2 Convert FASTA files to Nexus files
```
conda activate mtb_wgs_anno_env

CONSENSUS_FASTA_FILE=$WORKING_PATH/05_Concensus_Fasta
CONSENSUS_NEXUS_FILE=$WORKING_PATH/06_Concensus_Nexus

mkdir -p $CONSENSUS_NEXUS_FILE

cd $CONSENSUS_FASTA_FILE

for FASTA_FILE in $(find $CONSENSUS_FASTA_FILE -name "*.fasta" -exec basename "{}" \; | \
                    grep --invert-match 'merged_seq.fasta')
  do
    seqmagick convert \
      --input-format fasta \
      --output-format nexus \
      --alphabet dna \
      $FASTA_FILE $CONSENSUS_NEXUS_FILE/${FASTA_FILE/.fasta/.nex}
  done
```

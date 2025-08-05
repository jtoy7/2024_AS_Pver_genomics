# Pocillopora species assignment via UCE/exon mapping as in Oury et al. 2023

## Mapping
Will use UCE/Exon reference provided by Nico Oury. Contains 2068 consensus sequences that should work for mapping all Pocillopora species.
File is called `Pocillopora_UCEs_exons_2068_ref_sequences.fasta`.

### Prepare reference
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/UCE_exon_mapping

module load bwa
crun.bwa bwa index Pocillopora_UCEs_exons_2068_ref_sequences.fasta


module load samtools
crun.samtools samtools faidx Pocillopora_UCEs_exons_2068_ref_sequences.fasta


module load picard
crun.picard picard -Xmx4g CreateSequenceDictionary \
  REFERENCE=Pocillopora_UCEs_exons_2068_ref_sequences.fasta \
  OUTPUT=Pocillopora_UCEs_exons_2068_ref_sequences.dict
```


### Consolidate all trimmed fastq files into one location using symbolic links
Create symbolic links (shortcuts) for fastq files from all batches in directory `BASEDIR=/archive/barshis/barshislab/jtoy/`:
```bash
BASEDIR=/archive/barshis/barshislab/jtoy/

mkdir $BASEDIR/pver_gwas/UCE_exon_mapping/trimmed_fastq_links

# BATCH 1
FQDIR=/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/trimmed_fastq/
FQLIST=/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/sample_lists/fastq_trimmed_list_pver_pilot.txt
LINKDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/trimmed_fastq_links

for FILE in $(cat $FQLIST); do
  ln -s $FQDIR/$FILE $LINKDIR
done

## rename files to match batches 2 and 3
    ## syntax: rename <expression> <replacement> <file>
rename -v - _ *-*
rename -v _227H3WLT4_ _R24114_ *_227H3WLT4_*
rename -v _1_ _1_A_ *_1_*


# BATCH 2
FQDIR=$BASEDIR/pver_gwas/pver_gwas_batch2/trimmed_fastq
FQLIST=$BASEDIR/pver_gwas/pver_gwas_batch2/sample_lists/fastq_trimmed_list_pver_batch2.txt
LINKDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/trimmed_fastq_links

for FILE in $(cat $FQLIST); do
  ln -s $FQDIR/$FILE $LINKDIR
done


# BATCH 3
FQDIR=$BASEDIR/pver_gwas/pver_gwas_batch3/trimmed_fastq
FQLIST=$BASEDIR/pver_gwas/pver_gwas_batch3/sample_lists/fastq_trimmed_list_pver_batch3.txt
LINKDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/trimmed_fastq_links

for FILE in $(cat $FQLIST); do
  ln -s $FQDIR/$FILE $LINKDIR
done
```





### Run mapping array script:
`UCE_exon_mapping_array.slurm`
```bash
#!/bin/bash
#SBATCH --job-name=UCE_exon_mapping_array_2025-07-28
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-4694%24
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=10

# Load necessary modules
module load container_env bwa

# Define variables
BASEDIR=/archive/barshis/barshislab/jtoy/
FASTQDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/trimmed_fastq_links
OUTDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam
SAMPLELIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/fastq_trimmed_list_pver_gwas_all.txt
FASTQ_SUFFIX_1=_f_paired_trim.fastq.gz
FASTQ_SUFFIX_2=_r_paired_trim.fastq.gz
REFBASENAME=UCEexon2068ref

# Get SLURM task ID's sample file
echo $SLURM_JOB_ID
SAMPLEFILE=$(head -n $SLURM_ARRAY_TASK_ID $SAMPLELIST | tail -n 1)
echo $SAMPLEFILE

# Extract metadata for read group
SAMPLE_ID=$(echo $SAMPLEFILE | cut -f1-4 -d"_")
POP_ID=$(echo $SAMPLEFILE | cut -f2 -d"_")
SEQ_ID=$(echo $SAMPLEFILE | cut -f5 -d"_")
LANE_ID=$(echo $SAMPLEFILE | cut -f7,8 -d"_")
PREP_ID=$(echo $SAMPLEFILE | cut -f6 -d"_")
PU=$LANE_ID  # Platform unit
SAMPLE_UNIQ_ID=${SAMPLE_ID}_${SEQ_ID}_${PREP_ID}_${LANE_ID}
echo $SAMPLE_UNIQ_ID

# Set output file path
SAMPLEOUT=$OUTDIR/$SAMPLE_UNIQ_ID

# Construct read group string for bwa mem
RG="@RG\tID:${SAMPLE_UNIQ_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:${SAMPLE_ID}${SEQ_ID}${PREP_ID}\tPU:${PU}"

# Run bwa mem on fastq files directly
crun.bwa bwa mem -t 4 -R "${RG}" \
    $BASEDIR/pver_gwas/UCE_exon_mapping/references/Pocillopora_UCEs_exons_2068_ref_sequences.fasta \
    ${FASTQDIR}/${SAMPLE_UNIQ_ID}${FASTQ_SUFFIX_1} \
    ${FASTQDIR}/${SAMPLE_UNIQ_ID}${FASTQ_SUFFIX_2} \
    > ${SAMPLEOUT}_bwa_${REFBASENAME}.sam


## Change modules
module unload bwa
module load container_env gatk
GATK='crun.gatk gatk'
module load container_env samtools


## Query-sort for duplicate removal with GATK
# Run SortSam to sort by query name and convert to BAM
$GATK --java-options "-Xmx100G" SortSam \
  --INPUT ${SAMPLEOUT}_bwa_${REFBASENAME}.sam \
  --OUTPUT ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.bam \
  --SORT_ORDER queryname

# Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.bam \
  -O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.val \
  -M VERBOSE

rm ${SAMPLEOUT}_bwa_${REFBASENAME}.sam

## Mark and remove duplicates
$GATK --java-options "-Xmx100G" MarkDuplicates \
  -I ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.bam \
  -O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.bam \
  --METRICS_FILE ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dupstat.txt \
  --REMOVE_DUPLICATES true

## Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.bam \
  -O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.val \
  -M VERBOSE

#remove old BAM file
rm ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.bam

## Run SortSam to sort by coordinate for downstream processing
$GATK --java-options "-Xmx100G" SortSam \
  --INPUT ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.bam \
  --OUTPUT ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.bam \
  --SORT_ORDER coordinate

# Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.bam \
  -O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.val \
  -M VERBOSE

# Index final coordinate-sorted bam for use in GATK3 IndelRealigner
crun.samtools samtools index -@ 4 ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.bam

#remove old BAM file
rm ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.bam
```

Some jobs did not complete successfully because Java ran out of memory. I'm guessing this is because bwa is more memory intensive than bowtie2.

<br>

To get a list of which files need to be rerun with more memory, I ran:
```bash
cat 4630805*.out | grep -B 6 "There is insufficient memory" | grep "2024_" > ../sample_lists/mapping_rerun_list.txt
```
This lists **75 samples** that did not complete the job because of insufficeint memory.

Then removed the existing output files for these input files so they didn't interfere with the reruns:
```bash
for SAMPLE in `cat $BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/mapping_rerun_list.txt`; do
  echo $SAMPLE
  rm $BASEDIR/pver_gwas/UCE_exon_mapping/bam/${SAMPLE}'_bwa'*
done
```

I double checked the number of bams and there are actually **76 missing**.

So I compared a list of samples with the "insufficent memory" error in their .out files (`mapping_rerun_list.txt`), a list of samples with FAILED job status according to SLURM (`failed_prefixed.txt`), and a list of samples that are missing .bam files (`missing_prefixes.txt`):

```bash
diff <(sort missing_prefixes.txt) <(sort sample_lists/mapping_rerun_list.txt)
```
```
1d0
< 2024_ALOF_Pver_05_1_B_R24196_L005
```

<br>

```bash
diff <(sort missing_prefixes.txt) <(sort failed_prefixes.txt)
```
```
27d26
< 2024_ALOF_Pver_30_1_A_R24193_L006
```

<br>

```bash
diff <(sort sample_lists/mapping_rerun_list.txt) <(sort failed_prefixes.txt)
```
```
0a1
> 2024_ALOF_Pver_05_1_B_R24196_L005
26d26
< 2024_ALOF_Pver_30_1_A_R24193_L006
```

It turns out **job 93** (corresponding to the data file for `2024_ALOF_Pver_05_1_B_R24196_L005`) failed but for a reason other than insufficient memory. This was manually added to the list of samples to rerun.

Weirdly, **job 386** (corresponding to the data file for `2024_ALOF_Pver_30_1_A_R24193_L006`) was logged by SLURM as COMPLETED despite having an insufficient memory error in the .out file.

Then reran mapping array script for samples in `mapping_rerun_list.txt`, using greater memory allocation:
`UCE_exon_mapping_array_reun76samples.slurm`
```bash
#!/bin/bash
#SBATCH --job-name=UCE_exon_mapping_array_2025-08-04
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-76%12
#SBATCH --ntasks=1
#SBATCH --mem=190G
#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=20

# Load necessary modules
module load container_env bwa

# Define variables
BASEDIR=/archive/barshis/barshislab/jtoy/
FASTQDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/trimmed_fastq_links
OUTDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam
SAMPLELIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/mapping_rerun_list.txt
FASTQ_SUFFIX_1=_f_paired_trim.fastq.gz
FASTQ_SUFFIX_2=_r_paired_trim.fastq.gz
REFBASENAME=UCEexon2068ref

# Get SLURM task ID's sample file
echo $SLURM_JOB_ID
SAMPLEFILE=$(head -n $SLURM_ARRAY_TASK_ID $SAMPLELIST | tail -n 1)
echo $SAMPLEFILE

# Extract metadata for read group
SAMPLE_ID=$(echo $SAMPLEFILE | cut -f1-4 -d"_")
POP_ID=$(echo $SAMPLEFILE | cut -f2 -d"_")
SEQ_ID=$(echo $SAMPLEFILE | cut -f5 -d"_")
LANE_ID=$(echo $SAMPLEFILE | cut -f7,8 -d"_")
PREP_ID=$(echo $SAMPLEFILE | cut -f6 -d"_")
PU=$LANE_ID  # Platform unit
SAMPLE_UNIQ_ID=${SAMPLE_ID}_${SEQ_ID}_${PREP_ID}_${LANE_ID}
echo $SAMPLE_UNIQ_ID

# Set output file path
SAMPLEOUT=$OUTDIR/$SAMPLE_UNIQ_ID

# Construct read group string for bwa mem
RG="@RG\tID:${SAMPLE_UNIQ_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:${SAMPLE_ID}${SEQ_ID}${PREP_ID}\tPU:${PU}"

# Run bwa mem on fastq files directly
crun.bwa bwa mem -t 20 -R "${RG}" \
    $BASEDIR/pver_gwas/UCE_exon_mapping/references/Pocillopora_UCEs_exons_2068_ref_sequences.fasta \
    ${FASTQDIR}/${SAMPLE_UNIQ_ID}${FASTQ_SUFFIX_1} \
    ${FASTQDIR}/${SAMPLE_UNIQ_ID}${FASTQ_SUFFIX_2} \
    > ${SAMPLEOUT}_bwa_${REFBASENAME}.sam


## Change modules
module unload bwa
module load container_env gatk
GATK='crun.gatk gatk'
module load container_env samtools


## Query-sort for duplicate removal with GATK
# Run SortSam to sort by query name and convert to BAM
$GATK --java-options "-Xmx190G" SortSam \
  --INPUT ${SAMPLEOUT}_bwa_${REFBASENAME}.sam \
  --OUTPUT ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.bam \
  --SORT_ORDER queryname

# Run validation of BAM file
$GATK --java-options "-Xmx190G" ValidateSamFile \
  -I ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.bam \
  -O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.val \
  -M VERBOSE

rm ${SAMPLEOUT}_bwa_${REFBASENAME}.sam

## Mark and remove duplicates
$GATK --java-options "-Xmx190G" MarkDuplicates \
  -I ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.bam \
  -O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.bam \
  --METRICS_FILE ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dupstat.txt \
  --REMOVE_DUPLICATES true

## Run validation of BAM file
$GATK --java-options "-Xmx190G" ValidateSamFile \
  -I ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.bam \
  -O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.val \
  -M VERBOSE

#remove old BAM file
rm ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted.bam

## Run SortSam to sort by coordinate for downstream processing
$GATK --java-options "-Xmx190G" SortSam \
  --INPUT ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.bam \
  --OUTPUT ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.bam \
  --SORT_ORDER coordinate

# Run validation of BAM file
$GATK --java-options "-Xmx190G" ValidateSamFile \
  -I ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.bam \
  -O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.val \
  -M VERBOSE

# Index final coordinate-sorted bam for use in GATK3 IndelRealigner
crun.samtools samtools index -@ 20 ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.bam

#remove old BAM file
rm ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup.bam
```

<br>


### Summarize alignment metrics

Run CollectAlignmentSummaryMetrics from picard to get mapping summary for each file.
`CollectAlignmentSummaryMetrics_array.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=CollectAlignmentSummaryMetrics_array_2025-08-04
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-4694%42
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=6

# Load necessary modules
module load container_env gatk
GATK='crun.gatk gatk'

# Define variables
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam
SAMPLELIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/fastq_trimmed_list_pver_gwas_all.txt
REFERENCE=$BASEDIR/pver_gwas/UCE_exon_mapping/references/Pocillopora_UCEs_exons_2068_ref_sequences.fasta
REFBASENAME=UCEexon2068ref

# Get SLURM task ID's sample file
echo $SLURM_JOB_ID
SAMPLEFILE=$(head -n $SLURM_ARRAY_TASK_ID $SAMPLELIST | tail -n 1)
echo $SAMPLEFILE

# Extract metadata for read group
SAMPLE_ID=$(echo $SAMPLEFILE | cut -f1-4 -d"_")
POP_ID=$(echo $SAMPLEFILE | cut -f2 -d"_")
SEQ_ID=$(echo $SAMPLEFILE | cut -f5 -d"_")
LANE_ID=$(echo $SAMPLEFILE | cut -f7,8 -d"_")
PREP_ID=$(echo $SAMPLEFILE | cut -f6 -d"_")
PU=$LANE_ID  # Platform unit
SAMPLE_UNIQ_ID=${SAMPLE_ID}_${SEQ_ID}_${PREP_ID}_${LANE_ID}
echo $SAMPLE_UNIQ_ID

# Set output file path
SAMPLEOUT=$BAMDIR/$SAMPLE_UNIQ_ID

# Set input bam file name
BAMFILE=${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.bam


# Run CollectAlignmentSummaryMetrics
$GATK --java-options "-Xmx30G" CollectAlignmentSummaryMetrics \
-R $REFERENCE \
-I $BAMFILE \
-O ${SAMPLEOUT}_bwa_${REFBASENAME}_qsorted_dedup_coordsorted.bam_alignment_metrics.txt
```

Parse alignment metrics for each file into a summary table:
`parse_CollectAlignmentSummaryMetrics_output.sh`
```bash

```

Check duplicate rate of each library:
`summarize_dedup.sh`
```bash

```

Remove Ahya bams into a separate directory:
```bash
mkdir Ahya_bams
mv *Ahya*.bam Ahya_bams
```

### Merge bams from same extraction
Run script `merge_bams_by_sample_array.slurm` to merge bams from data batches 2 and 3 (4 + 8 = 12 lanes per library) by sample:
```bash
#!/bin/bash

#SBATCH --job-name=merge_bams_by_sample_2025-08-05
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-397%42
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=6

## Load modules
module load container_env samtools

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam
OUTDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/merged_bams
SAMPLELIST=($(ls $BAMDIR/*P*.bam | sed -E 's/.*\/2024_([A-Z]{3}._P[a-z]{3}_[0-9A-Za-z]+_[12])_.*/\1/' | sort | uniq))
REFBASENAME=UCEexon2068ref


# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"


# Get the sample name for current array task
SAMPLE="${SAMPLELIST[$SLURM_ARRAY_TASK_ID-1]}"

# Keep record of sample file
echo $SAMPLE

# List the BAM files for the current sample
BAMFILES=$(ls $BAMDIR/'2024_'$SAMPLE'_'*'.bam')

# Keep record of BAM files
echo $BAMFILES

# Keep a record of the Job ID
echo $SLURM_JOB_ID

# Define the output merged BAM file name
MERGEDBAM=$OUTDIR/'2024_'$SAMPLE'_'$REFBASENAME'_merged.bam'


# Merge the BAM files for the sample
echo "Merging BAM files for $SAMPLE into $MERGEDBAM..."
crun.samtools samtools merge -f -@ 6 "$MERGEDBAM" $BAMFILES


echo "Merging complete for $SAMPLE!"
```

```bash
sbatch merge_bams_by_sample_array.slurm
```

### Make list of merged bam files
```bash
cd $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams

ls *merged.bam > ../sample_lists/merged_bams_list.txt
```
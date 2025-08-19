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
#!/bin/bash

METRICS_DIR=/archive/barshis/barshislab/jtoy/pver_gwas/UCE_exon_mapping/bam
OUTPUT_FILE="$METRICS_DIR"/alignment_summary_metrics_combined.tsv

> "$OUTPUT_FILE"

header_written=0

for METRIC_FILE in `ls "$METRICS_DIR"/*_alignment_metrics.txt`; do
    SAMPLE_ID=$(basename "$METRIC_FILE" | cut -d'_' -f1-8)

    if [[ $header_written -eq 0 ]]; then
        HEADER_LINE=$(head -n 7 "$METRIC_FILE" | tail -n 1)
        echo -e "SAMPLE\t$HEADER_LINE" >> "$OUTPUT_FILE"
        header_written=1
    fi

    METRIC_LINES=$(head -n 10 "$METRIC_FILE" | tail -n 3)

    while IFS= read -r line; do
        echo -e "$SAMPLE_ID\t$line" >> "$OUTPUT_FILE"
    done <<< "$METRIC_LINES"
done

echo "Combined metrics saved to $OUTPUT_FILE"
```
The following files had no alignments, so they only output one CATEGORY = "UNPAIRED" line per file:
2024_AOAA_Pver_03_1_A_R24193_L005_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam_alignment_metrics.txt
2024_AOAA_Pver_03_1_A_R24196_L003_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam_alignment_metrics.txt
2024_AOAA_Pver_03_1_A_R24196_L004_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam_alignment_metrics.txt

This is why the `alignment_summary_metrics_combined.tsv` file is 6 lines shorter than expected.

<br>

Check duplicate rate of each library:
`summarize_dedup.sh`
```bash
#!/bin/bash

# Create output file
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam
OUTFILE=$BAMDIR/dupstat_summary.tsv


# Change working directory
cd $BAMDIR

# Print header
echo -e "FILE\tLIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE" > $OUTFILE

# Loop through all dupstat files
for FILE in `ls *dupstat.txt`; do
    # Extract the second line after the "## METRICS CLASS" line (i.e., the actual data)
    DATA=$(awk '/^## METRICS CLASS/ {getline; getline; print}' "$FILE")

    # Print data to outfile
    echo -e "$FILE\t$DATA" >> $OUTFILE

done
```

```bash
less -S dupstat_summary.tsv
```
```
FILE    LIBRARY READ_PAIRS_EXAMINED     PERCENT_DUPLICATION
2024_VATI_Pver_33_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_331A     1436337 0.744334
2024_VATI_Pver_36_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_361A     1581278 0.743764
2024_VATI_Pver_35_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_351A     1609489 0.70955
2024_VATI_Pver_30_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_301A     3409757 0.703349
2024_VATI_Pver_25_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_251A     3759393 0.692728
2024_VATI_Pver_22_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_221A     5224472 0.690568
2024_VATI_Pver_34_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_341A     5181677 0.688023
2024_VATI_Pver_24_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_241A     4882890 0.687029
2024_VATI_Pver_27_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_271A     3400378 0.686897
2024_VATI_Pver_26_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_261A     3157198 0.686232
2024_VATI_Pver_29_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_291A     5377729 0.683454
2024_VATI_Pver_32_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_321A     5650777 0.681824
2024_FTEL_Pver_14_1_B_R24196_L006_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_FTEL_Pver_141B     190639  0.666475
2024_VATI_Pver_28_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_281A     4158454 0.664436
2024_VATI_Pver_31_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_311A     5404326 0.653372
2024_VATI_Pver_23_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_231A     5567926 0.650706
2024_VATI_Pver_21_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_211A     6100682 0.649143
2024_OFU3_Pver_24_1_A_R24196_L002_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_241A     0       0.642857
2024_FTEL_Pver_14_1_B_R24196_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_FTEL_Pver_141B     189575  0.641471
2024_FTEL_Pver_14_1_B_R24196_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_FTEL_Pver_141B     189670  0.639482
...
2024_LEON_Ahya_21_1_B_R24196_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Ahya_211B     123537  0.235735
2024_LEON_Ahya_21_1_B_R24193_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Ahya_211B     121930  0.233313
2024_LEON_Ahya_21_1_B_R24196_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Ahya_211B     122547  0.233155
2024_LEON_Ahya_21_1_B_R24193_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Ahya_211B     120430  0.232021
2024_VATI_Pver_17_1_B_R24196_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_171B     355460  0.231093
2024_VATI_Pver_17_1_B_R24193_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_171B     347255  0.228347
2024_VATI_Pver_17_1_B_R24196_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_171B     354720  0.227863
2024_VATI_Pver_17_1_B_R24193_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_171B     353858  0.226652
2024_FMAL_Pspp_Extra1_1_A_R24196_L003_bwa_UCEexon2068ref_qsorted_dupstat.txt    2024_FMAL_Pspp_Extra11A 76      0.22488
2024_VATI_Pver_17_1_B_R24196_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_VATI_Pver_171B     348098  0.22313
2024_OFU3_Pver_12_1_A_R24196_L002_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_121A     85      0.221538
2024_OFU3_Pver_12_1_A_R24193_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_121A     114     0.21875
2024_FMAL_Pspp_Extra1_1_A_R24196_L004_bwa_UCEexon2068ref_qsorted_dupstat.txt    2024_FMAL_Pspp_Extra11A 82      0.218605
2024_FMAL_Pspp_Extra1_1_A_R24196_L001_bwa_UCEexon2068ref_qsorted_dupstat.txt    2024_FMAL_Pspp_Extra11A 74      0.207373
2024_OFU3_Pver_16_1_A_R24196_L001_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_161A     1       0.2
2024_FMAL_Pspp_Extra1_1_A_R24196_L002_bwa_UCEexon2068ref_qsorted_dupstat.txt    2024_FMAL_Pspp_Extra11A 79      0.186364
2024_LEON_Pver_03_1_B_R24196_L006_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_031B     97065   0.178101
2024_FMAL_Pspp_Extra1_1_A_R24193_L006_bwa_UCEexon2068ref_qsorted_dupstat.txt    2024_FMAL_Pspp_Extra11A 75      0.177419
2024_LEON_Pver_03_1_B_R24196_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_031B     97983   0.165869
2024_LEON_Pver_03_1_B_R24196_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_031B     98746   0.164037
2024_FMAL_Pspp_Extra1_1_A_R24193_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt    2024_FMAL_Pspp_Extra11A 66      0.163366
2024_LEON_Pver_03_1_B_R24193_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_031B     99020   0.161474
2024_LEON_Pver_03_1_B_R24196_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_031B     97216   0.160803
2024_LEON_Pver_03_1_B_R24193_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_031B     95774   0.158918
2024_OFU3_Pver_27_1_A_R24196_L001_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_271A     1       0.142857
2024_OFU3_Pver_24_1_A_R24196_L001_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_241A     1       0.142857
2024_OFU3_Pver_16_1_A_R24193_L006_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_161A     2       0.111111
2024_OFU3_Pver_16_1_A_R24196_L004_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_161A     2       0
2024_OFU3_Pver_16_1_A_R24196_L002_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_161A     2       0
2024_OFU3_Pver_14_1_A_R24196_L004_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_141A     0       0
2024_OFU3_Pver_14_1_A_R24196_L003_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_141A     1       0
2024_OFU3_Pver_14_1_A_R24196_L002_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_141A     0       0
2024_OFU3_Pver_14_1_A_R24196_L001_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_141A     0       0
2024_OFU3_Pver_14_1_A_R24193_L006_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_141A     0       0
2024_OFU3_Pver_14_1_A_R24193_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_OFU3_Pver_141A     0       0
2024_LEON_Pver_11_1_B_R24196_L006_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_111B     0       0
2024_LEON_Pver_11_1_B_R24196_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_111B     0       0
2024_LEON_Pver_11_1_B_R24193_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_111B     0       0
2024_LEON_Pver_11_1_B_R24193_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Pver_111B     0       0
2024_LEON_Ahya_21_1_A_R24196_L003_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Ahya_211A     1       0
2024_LEON_Ahya_21_1_A_R24193_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_LEON_Ahya_211A     0       0
2024_AOAA_Pver_03_1_B_R24196_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031B     0       0
2024_AOAA_Pver_03_1_B_R24196_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031B     0       0
2024_AOAA_Pver_03_1_B_R24196_L006_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031B     0       0
2024_AOAA_Pver_03_1_B_R24196_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031B     1       0
2024_AOAA_Pver_03_1_B_R24193_L008_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031B     0       0
2024_AOAA_Pver_03_1_B_R24193_L007_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031B     0       0
2024_AOAA_Pver_03_1_A_R24196_L004_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031A     0       0
2024_AOAA_Pver_03_1_A_R24196_L003_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031A     0       0
2024_AOAA_Pver_03_1_A_R24196_L002_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031A     0       0
2024_AOAA_Pver_03_1_A_R24196_L001_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031A     0       0
2024_AOAA_Pver_03_1_A_R24193_L006_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031A     0       0
2024_AOAA_Pver_03_1_A_R24193_L005_bwa_UCEexon2068ref_qsorted_dupstat.txt        2024_AOAA_Pver_031A     0       0
```
Most files have pretty high duplication rate around 40-60%, which is higher than for the hologenome mapping. I guess this is because there are fewer regions to map to with this reference, so reads can map to the same regions more easily?

<br>

Move Ahya bams into a separate directory:
```bash
mkdir Ahya_bams
mv *Ahya*.bam Ahya_bams
```
36 files moved.

<br>

Move Batch 1 bams into a separate directory (they have one file each so don't need to be merged):
```bash
mkdir batch1_bams
mv 2024_VATI_Pver_21_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_22_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_23_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_24_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_25_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_26_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_27_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_28_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_29_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_30_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_31_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_32_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_33_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_34_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_35_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
mv 2024_VATI_Pver_36_1_A_R24114_L008_bwa_UCEexon2068ref_qsorted_dedup_coordsorted.bam batch1_bams/
```
16 files moved.

<br>

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

### Dedup merged bams
`dedup_merged_bams_array.slurm`:
```bash
#!/bin/bash

#SBATCH --job-name dedup_merged_bams_array_2025-08-05
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-381%12
#SBATCH --ntasks=1
#SBATCH --mem=190G
#SBATCH --time 14-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env
module load container_env gatk
GATK='crun.gatk gatk'
module load samtools

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/merged_bams
OUTDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams
SAMPLELIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/merged_bams_list.txt

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Make directory
mkdir -p $OUTDIR


## Query-sort for duplicate removal with GATK
# Run SortSam to sort by query name and convert to BAM
$GATK --java-options "-Xmx190G" SortSam \
  --INPUT $BAMDIR/$SAMPLEFILE \
  --OUTPUT $BAMDIR/${SAMPLEFILE%.*}'_qsorted.bam' \
  --SORT_ORDER queryname

# Run validation of BAM file
$GATK --java-options "-Xmx190G" ValidateSamFile \
  -I $BAMDIR/${SAMPLEFILE%.*}'_qsorted.bam' \
  -O $BAMDIR/${SAMPLEFILE%.*}'_qsorted.val' \
  -M VERBOSE

## Mark and remove duplicates
$GATK --java-options "-Xmx190G" MarkDuplicates \
  -I $BAMDIR/${SAMPLEFILE%.*}'_qsorted.bam' \
  -O $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.bam' \
  --METRICS_FILE $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dupstat.txt' \
  --REMOVE_DUPLICATES true

## Run validation of BAM file
$GATK --java-options "-Xmx190G" ValidateSamFile \
  -I $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.bam' \
  -O $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.val' \
  -M VERBOSE

## Remove qsorted BAM files to make space
rm $BAMDIR/${SAMPLEFILE%.*}'_qsorted.bam'

## Run SortSam to sort by coordinate for downstream processing
$GATK --java-options "-Xmx190G" SortSam \
  --INPUT $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.bam' \
  --OUTPUT $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup_coordsorted.bam' \
  --SORT_ORDER coordinate

## Run validation of coordinate-sorted BAM file
$GATK --java-options "-Xmx190G" ValidateSamFile \
  -I $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup_coordsorted.bam' \
  -O $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup_coordsorted.val' \
  -M VERBOSE

# Index final coordinate-sorted bam for use in GATK3 IndelRealigner
crun.samtools samtools index -@ 20 $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup_coordsorted.bam'

## Remove qsorted deduped BAM files to make space
rm $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.bam'
```

<br>


Check duplication rate in merged bams:

`summarize_dedup_mergedbams.sh`
```bash
#!/bin/bash

# Create output file
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams
OUTFILE=$BAMDIR/dupstat_summary_mergedbams.tsv


# Change working directory
cd $BAMDIR

# Print header
echo -e "FILE\tLIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE" > $OUTFILE

# Loop through all dupstat files
for FILE in `ls *dupstat.txt`; do
    # Extract the second line after the "## METRICS CLASS" line (i.e., the actual data)
    DATA=$(awk '/^## METRICS CLASS/ {getline; getline; print}' "$FILE")

    # Print data to outfile
    echo -e "$FILE\t$DATA" >> $OUTFILE

done
```

```bash
 sort -k10,10 dupstat_summary_mergedbams.tsv | cut -f1,2,10 | less -S
```
```
FILE    LIBRARY PERCENT_DUPLICATION
2024_AOAA_Pver_03_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_AOAA_Pver_031A     0
2024_OFU3_Pver_14_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_141A     0
2024_OFU3_Pver_16_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_161A     0
2024_FMAL_Pspp_Extra1_1_UCEexon2068ref_merged_qsorted_dupstat.txt       2024_FMAL_Pspp_Extra11A 0.003988
2024_OFU3_Pver_12_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_121A     0.008197
2024_OFU3_Pver_18_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_181A     0.036847
2024_OFU3_Pver_27_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_271A     0.059701
2024_MALO_Pver_06_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_MALO_Pver_061A     0.096774
2024_FTEL_Pver_23_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_FTEL_Pver_231A     0.097307
2024_OFU3_Pver_24_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_241A     0.1
2024_MALO_Pver_10_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_MALO_Pver_101A     0.102068
2024_FALU_Pver_07_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_FALU_Pver_071A     0.102355
2024_FALU_Pver_05_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_FALU_Pver_051A     0.102444
2024_AOAA_Pver_19_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_AOAA_Pver_191A     0.103454
...
2024_OFU3_Pver_34_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_341A     0.156022
2024_FASA_Pver_16_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_FASA_Pver_161A     0.159285
2024_OFU3_Pver_33_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_331A     0.159395
2024_OFU6_Pver_25_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU6_Pver_251A     0.159652
2024_AOAA_Pver_21_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_AOAA_Pver_211A     0.161582
2024_MALO_Pver_28_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_MALO_Pver_281A     0.161631
2024_OFU3_Pver_02_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_021A     0.161936
2024_OFU3_Pver_20_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_201A     0.162667
2024_OFU6_Pver_31_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU6_Pver_311A     0.163593
2024_OFU6_Pver_10_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU6_Pver_101A     0.165318
2024_OFU3_Pver_07_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU3_Pver_071A     0.16784
2024_OFU6_Pver_09_2_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU6_Pver_092A     0.184142
2024_OFU6_Pver_22_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_OFU6_Pver_221A     0.184584
2024_LEON_Pver_07_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_LEON_Pver_071A     0.18908
2024_FASA_Pver_15_1_UCEexon2068ref_merged_qsorted_dupstat.txt   2024_FASA_Pver_151A     0.211009
```
Duplication rates were mostly in the 10-15% range.

<br>

### Realign reads around indels

Create list of batch1 (pilot) dedup bams:
```bash
cd $BASEDIR/pver_gwas/UCE_exon_mapping/bam/batch1_bams
ls *.bam > $BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/batch1_bams_list.txt
```

<br>

Filter bams to keep only aligned read pairs to reduce .bam size:
`filter_for_aligned_reads_array_batch1.slurm`
```bash
#!/bin/bash

#SBATCH --job-name filter_for_aligned_reads_array_batch1_2025-08-14
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 14-00:00:00
#SBATCH --cpus-per-task=16


## Load modules
module load container_env
module load container_env gatk
GATK='crun.gatk gatk'
module load samtools

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/batch1_bams
OUTDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/batch1_bams/filtered_bams
SAMPLELIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/batch1_bams_list.txt

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Make directory
mkdir -p $OUTDIR


# Filter BAM file to include only aligned read pairs
crun.samtools samtools view -@16 -b -f 3 $BAMDIR/$SAMPLEFILE > $OUTDIR/${SAMPLEFILE%.*}'_alignedpairs.bam'

## Run validation of filtered BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $OUTDIR/${SAMPLEFILE%.*}'_alignedpairs.bam' \
  -O $OUTDIR/${SAMPLEFILE%.*}'_alignedpairs.val' \
  -M VERBOSE


# Index filtered BAM for use in GATK3 IndelRealigner
crun.samtools samtools index -@ 16 $OUTDIR/${SAMPLEFILE%.*}'_alignedpairs.bam'
```

<br>

Create list of batch2/3 merged and deduped bams:
```bash
cd $BASEDIR/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams
ls *.bam > $BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/batch23_merged_deduped_bams_list.txt
```

<br>

`filter_for_aligned_reads_array_batch23.slurm`
```bash
#!/bin/bash

#SBATCH --job-name filter_for_aligned_reads_array_batch23_2025-08-14
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-381%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 14-00:00:00
#SBATCH --cpus-per-task=16


## Load modules
module load container_env
module load container_env gatk
GATK='crun.gatk gatk'
module load samtools

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams
OUTDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams/filtered_bams
SAMPLELIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/batch23_merged_deduped_bams_list.txt

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Make directory
mkdir -p $OUTDIR


# Filter BAM file to include only aligned read pairs
crun.samtools samtools view -@16 -b -f 3 $BAMDIR/$SAMPLEFILE > $OUTDIR/${SAMPLEFILE%.*}'_alignedpairs.bam'

## Run validation of filtered BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $OUTDIR/${SAMPLEFILE%.*}'_alignedpairs.bam' \
  -O $OUTDIR/${SAMPLEFILE%.*}'_alignedpairs.val' \
  -M VERBOSE


# Index filtered BAM for use in GATK3 IndelRealigner
crun.samtools samtools index -@ 16 $OUTDIR/${SAMPLEFILE%.*}'_alignedpairs.bam'
```

<br>

Make list of deduped merged bams (including the pilot samples):
```bash
ls /archive/barshis/barshislab/jtoy/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams/filtered_bams/*.bam > $BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/dedup_filtered_bam_fullpath.list
ls /archive/barshis/barshislab/jtoy/pver_gwas/UCE_exon_mapping/bam/batch1_bams/filtered_bams/*.bam >> $BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/dedup_filtered_bam_fullpath.list
```

<br>

Run RealignerTargetCreator and IndelRealigner:

`indel_realignment.slurm`
```bash
#!/bin/bash

#SBATCH --job-name indel_realignment_2025-08-14
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --time 14-00:00:00
#SBATCH --cpus-per-task=38

BASEDIR=/archive/barshis/barshislab/jtoy/
BAMLIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/dedup_filtered_bam_fullpath.list # Path to a list of merged, deduplicated, and overlap clipped bam files. Must be indexed. Full paths should be included. This file has to have a suffix of ".list".
REFERENCE=$BASEDIR/pver_gwas/UCE_exon_mapping/references/Pocillopora_UCEs_exons_2068_ref_sequences.fasta


## Realign around in-dels
# This is done across all samples at once!
module load container_env gatk/3.8

## Create list of potential in-dels
crun.gatk gatk3 -Xmx300g \
-T RealignerTargetCreator \
-R $REFERENCE \
-I $BAMLIST \
-o $BASEDIR'/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams/filtered_bams/target.intervals' \
-nt 38 # number of threads to use

## Run the indel realigner tool
crun.gatk gatk3 -Xmx300g \
-T IndelRealigner \
-R $REFERENCE \
-I $BAMLIST \
-targetIntervals $BASEDIR'/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams/filtered_bams/target.intervals' \
--nWayOut _realigned.bam
```
This ran successfully and RealignerTargetCreator identified 151 intervals for potential realignment, but because all samples were being processed simultaneously, there were 147 intervals that were skipped "because there were too many reads".

<br>

I decided to try rerunning it in a per-sample fasion, using an array job:

`indel_realignment_persample_array.slurm`
```bash
#!/bin/bash

#SBATCH --job-name indel_realignment_persample_array_2025-08-18
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-397%100
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time 14-00:00:00
#SBATCH --cpus-per-task=1

BASEDIR=/archive/barshis/barshislab/jtoy/
BAMLIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/dedup_filtered_bam_fullpath.list # Path to a list of merged, deduplicated, and overlap clipped bam files. Must be indexed. Full paths should be included. This file has to have a suffix of .list
REFERENCE=$BASEDIR/pver_gwas/UCE_exon_mapping/references/Pocillopora_UCEs_exons_2068_ref_sequences.fasta
OUTDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams/filtered_bams/realigned_bams

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $BAMLIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE


mkdir -p $OUTDIR


## Realign around in-dels
module load container_env
module load gatk/3.8


ulimit -a


## Run the indel realigner tool
crun.gatk gatk3 -Xmx16g \
-T IndelRealigner \
-R $REFERENCE \
-I $SAMPLEFILE \
-targetIntervals $BASEDIR'/pver_gwas/UCE_exon_mapping/bam/merged_bams/dedup_bams/filtered_bams/target.intervals' \
--maxReadsForRealignment 60000 \
--maxReadsInMemory 600000 \
-o $OUTDIR/$(basename ${SAMPLEFILE%.*})_realigned.bam

# index output
crun.gatk samtools index $OUTDIR/$(basename ${SAMPLEFILE%.*})_realigned.bam

# summarize alignment stats
crun.gatk samtools flagstat $OUTDIR/$(basename ${SAMPLEFILE%.*})_realigned.bam
```
This indeed resulted in no "too many reads" warnings.

<br>

Parse samtools flagstat output
```bash


```

## Call Variants with bcftools
`call_variants.slurm`
```bash
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMLIST=$BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/realigned_bam_fullpath.list
REFERENCE=$BASEDIR/pver_gwas/UCE_exon_mapping/references/Pocillopora_UCEs_exons_2068_ref_sequences.fasta
SITELIST=$BASEDIR/pver_gwas/UCE_exon_mapping/references/List_Phylo2023_17465SNP.txt
OUTDIR=$BASEDIR/pver_gwas/UCE_exon_mapping/vcf

module load bcftools

crun.bcftools bcftools mpileup -A -B -I -Ou -a AD,DP,SP,INFO/AD -f $REFERENCE \
	-R $SITELIST -b $BAMLIST \
  | crun.bcftools bcftools call -mv --threads 38 -Ob –o $OUTDIR/Calls_17465SNP_AS.bcf
```


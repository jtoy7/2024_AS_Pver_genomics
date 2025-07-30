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
#SBATCH --array=1-4694%64
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=4

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

Some jobs did not complete successfully because Java ran out of memory. I'm guessing this is because bwa is more memory intensive than bowtie2. To get a list of which files need to be rerun with more memory, I ran:
```bash
cat *.out | less | grep -B 6 "There is insufficient memory" | grep "2024_" | wc -l > ../sample_lists/mapping_rerun_list.txt
```

Then removed the existing output files for these input files so they didn't interfere with the reruns:
```bash
for SAMPLE in `cat $BASEDIR/pver_gwas/UCE_exon_mapping/sample_lists/mapping_rerun_list.txt`; do
  echo $SAMPLE
  rm $BASEDIR/pver_gwas/UCE_exon_mapping/bam/${SAMPLE}`_bwa`*
done
```

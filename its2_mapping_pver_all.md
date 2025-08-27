ITS2 mapping of 2024 P. verrucosa WGS data from American Samoa - all samples
================
**AUTHOR:** Jason A. Toy  
**DATE:** 2025-08-27 <br><br>

## Data Prep

Make symbolic links for all trimmed fastq files so that they are all in one place:
```bash
ln -s /cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/trimmed_fastq/*.fastq.gz /archive/barshis/barshislab/jtoy/pver_gwas/its2_mapping_pver_all/trimmed_fastq_links/
ln -s /archive/barshis/barshislab/jtoy/pver_gwas/pver_gwas_batch2/trimmed_fastq/*.fastq.gz /archive/barshis/barshislab/jtoy/pver_gwas/its2_mapping_pver_all/trimmed_fastq_links/
ln -s /archive/barshis/barshislab/jtoy/pver_gwas/pver_gwas_batch3/trimmed_fastq/*.fastq.gz /archive/barshis/barshislab/jtoy/pver_gwas/its2_mapping_pver_all/trimmed_fastq_links/
```

Make list of fasta file prefixes:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/its2_mapping_pver_all/trimmed_fastq_links/

ls *.fastq.gz | sed -E 's/_[fr]_paired_trim\.fastq\.gz$//' | sort | uniq > ../sample_lists/fastq_list_pver_its2.txt

wc -l ../sample_lists/fastq_list_pver_its2.txt
```
There are 4694 total pairs of fastq files.

<br>

Make combined sample table for all sequencing batches:
```bash

```

## Map to ITS2 Sequence Database
`its2_mapping_array.slurm`
```bash
#!/bin/bash

#SBATCH --job-name its2_mapping_array_pver_all_2025-08-27
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-4694%24
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --time 14-00:00:00
#SBATCH --cpus-per-task=10


## Load modules
module load container_env
module load bowtie2/2.4.1

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
FASTQDIR=$BASEDIR/pver_gwas/its2_mapping_pver_all/trimmed_fastq_links #path to trimmed fastq.gz files
OUTDIR=$BASEDIR/pver_gwas/its2_mapping_pver_all/its2_mapping
SAMPLELIST=$BASEDIR/pver_gwas/its2_mapping_pver_all/sample_lists/fastq_list_pver_its2.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_table_pver_its2.txt # Path to a sample table
FASTQ_SUFFIX_1=_f_paired_trim.fastq.gz # Suffix to trimmed fastq files. Use forward reads with paired-end data.
FASTQ_SUFFIX_2=_r_paired_trim.fastq.gz # Suffix to trimmed fastq files. Use reverse reads with paired-end data.
REFDIR=/cm/shared/courses/dbarshis/barshislab/jtoy/references/
REFBASENAME=Voolstra_SymPortal_published_div_20230315_Itrimmed

## Keep a record of the Job ID
echo $SLURM_JOB_ID
echo $SLURM_ARRAY_TASK_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library. This is for the naming of trimmed/processed files
SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 7`
PREP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 8`
LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 11`
SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$PREP_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 11` # Define platform unit (PU), which is the lane number


echo $SAMPLE_UNIQ_ID

## Create OUTDIR if it doesn't already exist
mkdir -p $OUTDIR

## Define the output path and file prefix
SAMPLEOUT=$OUTDIR/$SAMPLE_UNIQ_ID


## Run bowtie2
crun.bowtie2 bowtie2 -p 10 --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:${SAMPLE_ID}_${SEQ_ID}_${PREP_ID} --rg PU:$PU --rg PL:ILLUMINA \
        --local -x $REFDIR/$REFBASENAME -1 $FASTQDIR/$SAMPLE_UNIQ_ID$FASTQ_SUFFIX_1 -2 $FASTQDIR/$SAMPLE_UNIQ_ID$FASTQ_SUFFIX_2 \
        -S $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'

# Note that the -k argument is not used here, but rather the default behavior which searches for distinct, valid alignments for each read. When it finds a valid alignment, it continues looking for alignments that are nearly as good or better. The best alignment found is reported (randomly selected from among best if tied). bowtie2 also reports unaligned reads by default with the YT:Z:UP code, and we need to keep these to maintain proper mate-pair information. Otherwise the SAM/BAM files will not pass validation.


## Add version number to bowtie2 command header line in SAM file (This may not be necessary for your version of bowtie, but it was for mine). If your BAM validations below throw errors because of a missing version number in the BAM header, uncomment the line below:
sed -i 's/PN:bowtie2\tVN:/PN:bowtie2\tVN:2.4.1/' $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'
    

## Change modules
module unload bowtie2
module load container_env gatk
GATK='crun.gatk gatk'

## Query-sort for duplicate removal with GATK
# Run SortSam to sort by query name and convert to BAM
$GATK --java-options "-Xmx110G" SortSam \
  --INPUT $SAMPLEOUT'_bt2_'$REFBASENAME'.sam' \
  --OUTPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted.bam' \
  --SORT_ORDER queryname

# Run validation of BAM file
$GATK --java-options "-Xmx110G" ValidateSamFile \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted.val' \
  -M VERBOSE

#remove SAM file
rm $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'

## Mark and remove duplicates
$GATK --java-options "-Xmx110G" MarkDuplicates \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted_dedup.bam' \
  --METRICS_FILE $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted_dupstat.txt' \
  --REMOVE_DUPLICATES true

## Run validation of BAM file
$GATK --java-options "-Xmx110G" ValidateSamFile \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted_dedup.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted_dedup.val' \
  -M VERBOSE

## If each sample has only one BAM file, you will normally want to convert back to coordinate sorted at this point, but because in this case I'm first going to merge multiple bams from the same sample, it is unnecessary at this point.

## Index BAM file
module load samtools
crun.samtools samtools index -@ 10 $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted_dedup.bam'
```
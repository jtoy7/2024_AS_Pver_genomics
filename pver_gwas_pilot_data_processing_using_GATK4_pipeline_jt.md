Processing of 2024 P. verrucosa WGS data from American Samoa using GATK4 tools and adapted best practices
================
Jason Toy
2024-08-08

## 0. A note on GATK syntax

GATK4 and Picard (both created and maintained by the Broad Institute)
are installed as java packages (GATK4 now comes with Picard included),
but can now be invoked without the java -jar command. Instead they are
invoked with the simple “gatk” command and any java arguments can
immediately follow the command, placed in quotation marks. The general
syntax looks like this:

``` bash
gatk [--java-options "jvm args like -Xmx4G go here"] ToolName [GATK args go here]
```

 

An example command:

``` bash
gatk --java-options "-Xmx100G" HaplotypeCaller -R reference.fasta -I input.bam -O output.vcf
```

Note that this is different from traditional Picard syntax which instead
uses the `Input=input.bam` syntax.  

## 1. Start with raw mapped SAMs from Bowtie2

Do **not** sort output SAMs with samtools! (see Step 2)  

Run ValidateSamFile before continuing to the next step:

``` bash
crun.gatk gatk ValidateSamFile \
  -I input.bam \
  -MODE SUMMARY
```

 

Array script for validation:

``` bash
#!/bin/bash
#SBATCH --job-name validate_bams_array_2024-08-09
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/reheadered_bams_list.txt
GATK='crun.gatk gatk'


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $BASEDIR'/pver_gwas_pilot/bam/'$SAMPLEBAM \
  -O $BASEDIR'/pver_gwas_pilot/bam/'${SAMPLEBAM%.*}'.val' \
  -M VERBOSE

echo 'done-zo woot!'
```

 

To obtain more detailed information about specific warnings and errors,
rerun the command with the MODE option set to “VERBOSE”. See
documentation for ValidateSamFile
[here](https://gatk.broadinstitute.org/hc/en-us/articles/360036351112-ValidateSamFile-Picard).

It is a recommended best practice to run ValidateSamFile on alignment
files after every step of processes, so that errors can be identified
sooner rather than later.  

## 2. Sort (by query) and convert to BAM using SortSam (Picard/GATK4)

Picard’s MarkDuplicates can more exhaustively look for duplicates if the
file is sorted by read-name (query-sorted). Query-sorted BAM/SAM files
are sorted based upon their read names and ordered lexiographically.
Coordinate-sorted BAM/SAM files are sorted by their aligned sequence
name (chromosome/linkage group/scaffold) and position. Picard can mark
and remove duplicates in either coordinate-sorted or query-sorted
BAM/SAM files, however, if the alignments are query-sorted it can test
secondary alignments for duplicates and mark unmapped mates of reads
identified as duplicates
([source](https://hbctraining.github.io/variant_analysis/lessons/04_alignment_file_processing.html)).
From Picard’s [MarkDuplicates
documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard):

“The program can take either coordinate-sorted or query-sorted inputs,
however the behavior is slightly different. When the input is
coordinate-sorted, unmapped mates of mapped records and
supplementary/secondary alignments are not marked as duplicates.
However, when the input is query-sorted (actually query-grouped), then
unmapped mates and secondary/supplementary reads are not excluded from
the duplication test and can be marked as duplicate reads.”

From experience, it seems like running MarkDuplicates on
coordinate-sorted BAMs (done with Samtools sort) can results in
“ERROR:MATE_NOT_FOUND” errors when running ValidateSamFile to check file
formatting and integrity. Based on the above information, it seems
likely that unmapped mates of mapped duplicates that are not filtered
out are the cause. Because of this, it makes sense to first query-sort
our SAM file and convert it to a BAM file with SamSort if the data will
eventually be run through GATK.  

``` bash
crun.gatk gatk SortSam \
  --INPUT $SAM_FILE \
  --OUTPUT $QUERY_SORTED_BAM_FILE \
  --SORT_ORDER queryname
```

 

Array script for query-sorting:

``` bash
#!/bin/bash
#SBATCH --job-name query_sort_array_2024-08-09
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/reheadered_bams_list.txt
GATK='crun.gatk gatk'


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run SortSam to sort by query name
$GATK --java-options "-Xmx100G" SortSam \
  --INPUT $BASEDIR'/pver_gwas_pilot/bam/'$SAMPLEBAM \
  --OUTPUT $BASEDIR'/pver_gwas_pilot/bam/'${SAMPLEBAM%.*}'_qsorted.bam' \
  --SORT_ORDER queryname

echo 'done-zo woot!'
```

### Validate BAMS before moving to next step:

Array script for validations:

``` bash
#!/bin/bash
#SBATCH --job-name validate_bams_array_2024-08-09
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/reheadered_qsorted_bams_list.txt
GATK='crun.gatk gatk'


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $BASEDIR'/pver_gwas_pilot/bam/'$SAMPLEBAM \
  -O $BASEDIR'/pver_gwas_pilot/bam/'${SAMPLEBAM%.*}'.val' \
  -M VERBOSE

echo 'done-zo woot!'
```

 

## 3. Mark and remove duplicates using MarkDuplicates

``` bash
crun.gatk gatk MarkDuplicates \
  --INPUT $QUERY_SORTED_BAM_FILE \
  --OUTPUT $DEDUPED_BAM_FILE \
  --METRICS_FILE $METRICS_FILE \
  --REMOVE_DUPLICATES true
```

 

Array script for marking and removing duplicates

``` bash
#!/bin/bash
#SBATCH --job-name mark_remove_duplicates_array_2024-08-10
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/reheadered_qsorted_bams_list.txt
GATK='crun.gatk gatk'


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Mark and remove duplicates
$GATK --java-options "-Xmx100G" MarkDuplicates \
  -I $BASEDIR'/pver_gwas_pilot/bam/'$SAMPLEBAM \
  -O $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/'${SAMPLEBAM%.*}'_dedup.bam' \
  --METRICS_FILE $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/'${SAMPLEBAM%.*}'_dupstat.txt' \
  --REMOVE_DUPLICATES true

echo 'done-zo woot!'
```

 

### Validate BAMS before moving to next step:

Array script for validations:

``` bash
#!/bin/bash
#SBATCH --job-name validate_dedup_bams_array_2024-08-10
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_bams_list.txt
GATK='crun.gatk gatk'


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/'$SAMPLEBAM \
  -O $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/'${SAMPLEBAM%.*}'.val' \
  -M VERBOSE

echo 'done-zo woot!'
```

 

## 4. Coordinate-sort and index deduped BAM files

Downstream analyses will require BAMs to be coordinate-sorted (not
query-sorted) and indexed. Both can be done simultaneously with SortSam:

``` bash
crun.gatk gatk SortSam \
  --INPUT $DEDUPED_BAM_FILE \
  --OUTPUT $COORDINATE_SORTED_BAM_FILE \
  --SORT_ORDER coordinate \
  --CREATE_INDEX true
```

 

Array script for coordinate-sorting:

``` bash
#!/bin/bash
#SBATCH --job-name coordinate_sort_array_2024-08-12
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_bams_list.txt
GATK='crun.gatk gatk'


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run SortSam to sort by coordinate
$GATK --java-options "-Xmx100G" SortSam \
  --INPUT $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/'$SAMPLEBAM \
  --OUTPUT $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/'${SAMPLEBAM%.*}'_coordsorted.bam' \
  --SORT_ORDER coordinate
  --CREATE_INDEX true

echo 'done-zo woot!'
```

 

### Count alignments remaining post-dedup

``` bash
#!/bin/bash
#SBATCH --job-name count_postdedup_reads_array_2024-08-14
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env samtools

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_bams_coordsorted_list.txt

## Change working directory
cd $BASEDIR/pver_gwas_pilot/bam/dedup_bams2

## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


COUNT=`crun.samtools samtools view $SAMPLEBAM | wc -l`

echo -e "$SAMPLEBAM\t$COUNT" >> postdedup_read_counts.txt


echo "done-zo woot!"
```

 

## 5. Calculate depth/coverage statistics with CollectWgsMetrics:

``` bash
crun.gatk gatk CollectWgsMetrics \
    --INPUT $COORDINATE_SORTED_BAM_FILE \
    --OUTPUT $METRICS_OUTPUT_FILE \
    --REFERENCE_SEQUENCE $REFERENCE
```

 

Array script for CollectWgsMetrics:

``` bash
#!/bin/bash
#SBATCH --job-name CollectWgsMetrics_array_2024-08-12
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_bams_coordsorted_list.txt
GATK='crun.gatk gatk'
REFERENCE=$BASEDIR/references/genomes/combined_pver_cd_hologenome.fa


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run CollectWgsMetrics
$GATK --java-options "-Xmx100G" CollectWgsMetrics \
  --INPUT $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/'$SAMPLEBAM \
  --OUTPUT $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/'${SAMPLEBAM%.*}'_metrics.txt' \
  --REFERENCE_SEQUENCE $REFERENCE

echo 'done-zo woot!'
```

 

## 6. Extract only host (Pver) alignments

Filter bams using `select_host_alignments_array.slurm` script:

``` bash
#!/bin/bash

#SBATCH --job-name select_host_alignments_array_2024-08-12
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20



## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
OUTDIR=$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_bams_coordsorted_list.txt # Path to a bam list
REFBASENAME=PverCD
SCAFLIST=$BASEDIR/references/genomes/Pver_scaffold_names_singleline.txt


# Make new directory
mkdir $OUTDIR

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Define the output file name
HOSTOUT=${SAMPLEFILE%combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam}$REFBASENAME'_dedup_primary_minq20_mlen20_pver.bam'

# Load module
module load container_env samtools


#index sorted BAM
crun.samtools samtools index -@ 20 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/$SAMPLEFILE

# Extract mappings that are primary alignments only (no unmapped or secondary/supplementary reads), with mapping score > 20, mapping length (CIGAR) > 20, and only on host scaffolds. Then extract only reads that aligned concordantly (-f 2)
crun.samtools samtools view -b -F 260 -q 20 -m 20 -@ 28 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/$SAMPLEFILE `cat $SCAFLIST` | crun.samtools samtools view -f 2 -b -o $OUTDIR/$HOSTOUT

# Remove Cgoreaui and Dtrenchii sequence header lines, then reheader bam
crun.samtools samtools view -H $OUTDIR/$HOSTOUT | sed -e '/Cgoreaui/d' -e '/Dtrenchii/d' > $OUTDIR/'header_'$SLURM_ARRAY_TASK_ID'.sam'

crun.samtools samtools reheader $OUTDIR/'header_'$SLURM_ARRAY_TASK_ID'.sam' $OUTDIR/$HOSTOUT > $OUTDIR/${HOSTOUT%.*}'_reheadered.bam'

# Remove original bam file
rm $OUTDIR/$HOSTOUT

    # edited this to remove header lines for sym contigs and reran on 2025-03-31
```

 

### Count host reads

``` bash
#!/bin/bash
#SBATCH --job-name count_host_reads_array_2024-08-14
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env samtools

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/pver_bams_list.txt

## Change working directory
cd $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams

## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


COUNT=`crun.samtools samtools view $SAMPLEBAM | wc -l`

echo -e "$SAMPLEBAM\t$COUNT" >> pver_read_counts.txt


echo "done-zo woot!"
```

### Calculate percent host reads

```bash
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy

# first sort counts files
sort -k1,1 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/postdedup_read_counts.txt >  $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/postdedup_read_counts_sorted.txt
sort -k1,1 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_counts.txt > $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_counts_sorted.txt

# paste columns together
paste <(cut -f1 "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_counts_sorted.txt") \
      <(cut -f2 "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_counts_sorted.txt") \
      <(cut -f2 "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/postdedup_read_counts_sorted.txt") \
      > "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_and_total_read_counts.txt"

# calculate percentages
awk '{print $0, $2/$3}' \
    "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_and_total_read_counts.txt" \
    > "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_percents.txt"
```

| Sample                                                                                       | Post-processing_Pver_read_count_condordant_only | Post-processing_total_read_count | Percent_Pver_reads |
|----------------------------------------------------------------------------------------------|-------------------------------------------------|----------------------------------|--------------------|
| 2024_VATI_Pver_21_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 176042818                                       | 391304080                        | 0.449888           |
| 2024_VATI_Pver_22_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 168635188                                       | 382997872                        | 0.440303           |
| 2024_VATI_Pver_23_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 151224676                                       | 342971808                        | 0.440925           |
| 2024_VATI_Pver_24_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 137435262                                       | 320073616                        | 0.429386           |
| 2024_VATI_Pver_25_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 116004972                                       | 271652826                        | 0.427034           |
| 2024_VATI_Pver_26_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 110263344                                       | 253683358                        | 0.434649           |
| 2024_VATI_Pver_27_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 112565296                                       | 261762888                        | 0.430028           |
| 2024_VATI_Pver_28_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 116284158                                       | 276011764                        | 0.421301           |
| 2024_VATI_Pver_29_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 138842314                                       | 325343822                        | 0.426756           |
| 2024_VATI_Pver_30_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 99908544                                        | 231446692                        | 0.43167            |
| 2024_VATI_Pver_31_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 127183126                                       | 299238960                        | 0.425022           |
| 2024_VATI_Pver_32_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 149764782                                       | 351492214                        | 0.426083           |
| 2024_VATI_Pver_33_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 73224156                                        | 168179902                        | 0.435392           |
| 2024_VATI_Pver_34_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 134343792                                       | 315951282                        | 0.425204           |
| 2024_VATI_Pver_35_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 69922392                                        | 160741570                        | 0.434999           |
| 2024_VATI_Pver_36_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 72196466                                        | 165617090                        | 0.435924           |



## 7. Calculate depth/coverage statistics for host with CollectWgsMetricsWithNonZeroCoverage:

``` bash
crun.gatk gatk CollectWgsMetricsWithNonZeroCoverage \
    --INPUT $COORDINATE_SORTED_BAM_FILE \
    --OUTPUT $METRICS_OUTPUT_FILE \
    --CHART $CHART_OUTPUT_FILE \
    --REFERENCE_SEQUENCE $REFERENCE
```

 

Array script for CollectWgsMetricsWithNonZeroCoverage:

``` bash
#!/bin/bash
#SBATCH --job-name CollectWgsMetrics_pver_only_nonzero_array_2025-03-31
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/pver_bams_list.txt
GATK='crun.gatk gatk'
REFERENCE=$BASEDIR/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run CollectWgsMetrics
$GATK --java-options "-Xmx100G" CollectWgsMetricsWithNonZeroCoverage \
  --INPUT $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/pver_bams/'$SAMPLEBAM \
  --OUTPUT $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/pver_bams/'${SAMPLEBAM%.*}'_metrics_nonzero.txt' \
  --CHART $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/pver_bams/'${SAMPLEBAM%.*}'_metrics_nonzero_chart.pdf' \
  --REFERENCE_SEQUENCE $REFERENCE

echo 'done-zo woot!'
```

 

### Plot histograms of coverage depth and summarise coverage statistics using R script `summarize_coverage_stats.R`

First extract histogram data from output files and save to new files:

``` bash
for FILE in `ls *nonzero.txt`; do
    tail -n +12 $FILE > ${FILE%.*}_histdata.tsv
done
```


 

Then run R script:

``` r
setwd("/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/bam/dedup_bams2/pver_bams/")

file_list <- list.files(pattern = "histdata\\.tsv$")

for (FILE in file_list) {
  
  # extract sample name
  sample <- substr(FILE, 1, 34)
  print(sample)
  
  # load data
  histdat <- read_tsv(file = FILE, col_names = TRUE) %>% 
    mutate(coverage_count = coverage * count_NON_ZERO_REGIONS)

  # plot non-zero depth histogram
  cov_plot <- ggplot(histdat) + 
    geom_col(aes(x = coverage, y = count_NON_ZERO_REGIONS))

  ggsave(cov_plot, filename = paste0(sample, "_coverage_hist.pdf"), device = "pdf")
  
  # calculate mean depth
  mean_depth = sum(histdat$coverage_count)/sum(histdat$count_NON_ZERO_REGIONS)

  # calculate median depth  
  expanded_dataset <- rep(histdat$coverage, histdat$count_NON_ZERO_REGIONS)
  median_depth <- median(expanded_dataset)

  # calculate % of genome covered
  pver_genome_size <- 353416094
  coverage <- sum(histdat$count_NON_ZERO_REGIONS)/pver_genome_size

  # write stats to file
  header <- "Sample\tMean_Depth\tMedian_Depth\tProportion_Genome_Covered"
  
  if (!file.exists("coverage_summary.tsv") || file.info("coverage_summary.tsv")$size == 0) {
    # If the file doesn't exist or is empty, write the header and the first line of data
    write(header, file="coverage_summary.tsv")
    write(paste(sample, mean_depth, median_depth, coverage, sep="\t"), file="coverage_summary.tsv", append=TRUE)
  } else {
    # If the file exists and is not empty, append the new line of data
    write(paste(sample, mean_depth, median_depth, coverage, sep="\t"), file="coverage_summary.tsv", append=TRUE)
  }
}
```

 

Example depth histogram:  

VATI_PVER_21
![](./2024_VATI_Pver_21_1_227H3WLT4-L008_coverage_hist.png)

Summary table:

| Sample                             | Mean_Depth       | Median_Depth | Proportion_Genome_Covered |
|------------------------------------|------------------|--------------|---------------------------|
| 2024_VATI_Pver_21_1_227H3WLT4-L008 | 89.0063505367799 | 103          | 0.726539898321665         |
| 2024_VATI_Pver_22_1_227H3WLT4-L008 | 87.1223636380977 | 100          | 0.726553180116353         |
| 2024_VATI_Pver_23_1_227H3WLT4-L008 | 76.0735135633069 | 87           | 0.724353509492412         |
| 2024_VATI_Pver_24_1_227H3WLT4-L008 | 70.8338457798976 | 80           | 0.720388675904499         |
| 2024_VATI_Pver_25_1_227H3WLT4-L008 | 59.9018946033327 | 68           | 0.721541133324845         |
| 2024_VATI_Pver_26_1_227H3WLT4-L008 | 58.2398005777285 | 65           | 0.721126732276092         |
| 2024_VATI_Pver_27_1_227H3WLT4-L008 | 59.3436637907855 | 66           | 0.718120595832288         |
| 2024_VATI_Pver_28_1_227H3WLT4-L008 | 60.1630501174668 | 67           | 0.714952104020481         |
| 2024_VATI_Pver_29_1_227H3WLT4-L008 | 70.965564059926  | 81           | 0.719696955849441         |
| 2024_VATI_Pver_30_1_227H3WLT4-L008 | 52.252016861391  | 59           | 0.717566303021843         |
| 2024_VATI_Pver_31_1_227H3WLT4-L008 | 62.5461639984566 | 71           | 0.720830548820451         |
| 2024_VATI_Pver_32_1_227H3WLT4-L008 | 76.0627980932915 | 87           | 0.721541195574415         |
| 2024_VATI_Pver_33_1_227H3WLT4-L008 | 40.4523427718123 | 44           | 0.711358719277793         |
| 2024_VATI_Pver_34_1_227H3WLT4-L008 | 68.7201314517914 | 78           | 0.71907164759735          |
| 2024_VATI_Pver_35_1_227H3WLT4-L008 | 38.0867020241067 | 42           | 0.710261361215769         |
| 2024_VATI_Pver_36_1_227H3WLT4-L008 | 39.436154917118  | 43           | 0.7122055086716           |

 

## 8. **OPTIONAL** - Base quality score recalibration (BQSR):  (**SKIP FOR NOW**)

This step is performed per-sample and consists of applying machine
learning to detect and correct for patterns of systematic errors in the
base quality scores, which are confidence scores emitted by the
sequencer for each base. In a nutshell, it is a data pre-processing step
that detects systematic errors made by the sequencing machine when it
estimates the accuracy of each base call.

This step requires a database of known variants to be included as input,
which often does not exist for non-model species. In these cases it is
recommended to bootstrap your own set of variants by doing an initial
round of variant calling and using only the high-confidence variants as
the known set.

From the GATK Documentation:

“Here’s how you would bootstrap a set of known variants:

- First do an initial round of variant calling on your original,
  unrecalibrated data.
- Then take the variants that you have the highest confidence in and use
  that set as the database of known variants by feeding it as a VCF file
  to the BaseRecalibrator.
- Finally, do a real round of variant calling with the recalibrated
  data. These steps could be repeated several times until convergence.”



## 9. Call variants per sample with HaplotypeCaller (in GVCF mode)


First index the P. verrucosa reference fasta:

    module load samtools

    cd /cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1
    
    crun.samtools samtools faidx GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta


    module load gatk
    GATK='crun.gatk gatk'
    $GATK --java-options "-Xmx100G" CreateSequenceDictionary --REFERENCE GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta
         


Basic usage:
``` bash
crun.gatk gatk HaplotypeCaller \
    --INPUT $COORDINATE_SORTED_BAM_FILE \
    --OUTPUT $GVCF_OUTPUT_FILE \
    --REFERENCE_SEQUENCE $REFERENCE
    -ERC GVCF

# -ERC GVCF (--emit-ref-confidence)    Controls how variant calls and reference confidence values are emitted.
#                                      Ensures that the output is a gVCF (Genomic VCF), which includes reference blocks for non-variant sites.
```


Array script for HaplotypeCaller:

``` bash
#!/bin/bash
#SBATCH --job-name HaplotypeCaller_array_pver_2025-03-26
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20


## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/pver_bams_list.txt
GATK='crun.gatk gatk'
REFERENCE=$BASEDIR/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF_036669915.1_ASM3666991v2_genomic.fna


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Index BAM
$GATK --java-options "-Xmx100G" BuildBamIndex \
      --INPUT $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/pver_bams/'$SAMPLEBAM \
      --OUTPUT $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/pver_bams/'${SAMPLEBAM%.*}'.bai'


## Run HaplotypeCaller
$GATK --java-options "-Xmx100G" HaplotypeCaller \
  -I $BASEDIR'/pver_gwas_pilot/bam/dedup_bams2/pver_bams/'$SAMPLEBAM \
  -O $BASEDIR'/pver_gwas_pilot/gvcfs/'${SAMPLEBAM%.*}'.g.vcf.gz' \
  -R $REFERENCE \
  -ERC GVCF
  --native-pair-hmm-threads 28

echo 'done-zo woot!'
```


This version of the script can be used to check to see if .bai index files already exist before running BuildBamIndex:

```bash
#!/bin/bash
#SBATCH --job-name HaplotypeCaller_array_pver_2025-03-26
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=20

## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/pver_bams_list.txt
GATK='crun.gatk gatk'
REFERENCE=$BASEDIR/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta

## Get sample BAM filename
SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BAMLIST)
BAMFILE=$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/$SAMPLEBAM
BAIFILE=${BAMFILE%.*}.bai

echo "Slurm array task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample BAM: $SAMPLEBAM"

## Check if BAM index exists, and create it if missing
if [[ ! -f "$BAIFILE" ]]; then
  echo "Index file $BAIFILE not found. Creating index..."
  cd $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/
  $GATK --java-options "-Xmx100G" BuildBamIndex --INPUT "$BAMFILE"
else
  echo "Index file $BAIFILE already exists. Skipping indexing."
fi

## Run HaplotypeCaller
$GATK --java-options "-Xmx100G" HaplotypeCaller \
  -I $BAMFILE \
  -O $BASEDIR'/pver_gwas_pilot/gvcfs/'${SAMPLEBAM%.*}'.g.vcf.gz' \
  -R $REFERENCE \
  -ERC GVCF \
  --native-pair-hmm-threads 28

echo "done-zo woot!"
```
Run time: 21 hrs for the largest file

## 10. Consolidate per sample GVCFs with GenomicsDBImport

First need to create a tab-delimited map file with sample_name in the first column and /path/to/vcf in the second column:


Basic usage:
```bash
    gatk --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path my_database \
       --batch-size 50 \      # use if processing many samples simultaneously. Limits memory usage.
       -L chr1:1000-10000 \
       --sample-name-map cohort.sample_map \
       --tmp-dir /path/to/large/tmp \    # can be used to specify an alternate temporary storage location with sufficient space if needed
       --reader-threads 5
```

Slurm script:
```bash
#!/bin/bash

#SBATCH --job-name GenomicsDBImport_pver_2025-04-03
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=38

## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
GATK='crun.gatk gatk'

$GATK --java-options "-Xmx100g -Xms10g" \
   GenomicsDBImport \
   --genomicsdb-workspace-path $BASEDIR/pver_gwas_pilot/genomicsdb/ \
   --sample-name-map $BASEDIR/pver_gwas_pilot/gvcfs/pver_gwas_pilot_gvcf.sample_map \
   -L $BASEDIR/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/genome_regions.list \
```
Run time: 10.5 hrs for 16 samples. The tools imports by scaffold. There are 52 and the largest scaffolds (chromosomes) took about 1 hr to import. Processing time for non-chromosome scaffolds was negligible.

## 11. Joint genotyping with GenotypeGVCFs

This step is parallelized by genomic region to speed up processing time, since GenotypeGVCFs doesn't have internal parallelization options.

 `GenotypeGVCFs_array.slurm`
 ```bash
#!/bin/bash

#SBATCH --job-name GenotypeGVCFs_pver_2025-04-07
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-52%52
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=10

## Load modules
module load container_env gatk

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
REFERENCE=$BASEDIR/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta
OUTDIR=$BASEDIR/pver_gwas_pilot/vcf
GENDB=$BASEDIR/pver_gwas_pilot/genomicsdb
SCAFLIST=$BASEDIR/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/genome_regions.list
GATK='crun.gatk gatk'


SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SCAFLIST)


echo "Processing scaffold: $SCAF"

$GATK --java-options "-Xmx100g" GenotypeGVCFs \
   -R $REFERENCE \
   -V gendb://$GENDB \
   -L $SCAF \
   -O $OUTDIR/$SCAF_'genotypes.vcf.gz' \
```

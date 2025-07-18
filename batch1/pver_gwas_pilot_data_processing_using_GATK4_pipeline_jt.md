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
This actually counts all reads, not just mapped reads so I reran it as the version below:
```bash
#!/bin/bash
#SBATCH --job-name count_postdedup_mapped_reads_array_2025-06-30
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=16


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


COUNT=`crun.samtools samtools view -@16 -F 4 $SAMPLEBAM | wc -l`

echo -e "$SAMPLEBAM\t$COUNT" >> postdedup_mapped_read_counts.txt


echo "done-zo woot!"
```
Sort counts file:
```bash
sort -k1,1 postdedup_mapped_read_counts.txt > postdedup_mapped_read_counts_sorted.txt
```

```
2024_VATI_Pver_21_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      336872789
2024_VATI_Pver_22_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      328158377
2024_VATI_Pver_23_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      292872027
2024_VATI_Pver_24_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      271939942
2024_VATI_Pver_25_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      229791777
2024_VATI_Pver_26_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      216203805
2024_VATI_Pver_27_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      222965633
2024_VATI_Pver_28_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      235054474
2024_VATI_Pver_29_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      272361527
2024_VATI_Pver_30_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      195771705
2024_VATI_Pver_31_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      252517657
2024_VATI_Pver_32_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      299038006
2024_VATI_Pver_33_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      144087190
2024_VATI_Pver_34_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      268303671
2024_VATI_Pver_35_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      137841677
2024_VATI_Pver_36_1_227H3WLT4-L008_bt2_combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam      141689136
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
Percent of total reads:
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

<br>

Percent of total **mapped** reads:
```bash
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy

# first sort counts files
sort -k1,1 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/postdedup_mapped_read_counts.txt >  $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/postdedup_mapped_read_counts_sorted.txt
sort -k1,1 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_counts.txt > $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_counts_sorted.txt

# paste columns together
paste <(cut -f1 "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_counts_sorted.txt") \
      <(cut -f2 "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_read_counts_sorted.txt") \
      <(cut -f2 "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/postdedup_mapped_read_counts_sorted.txt") \
      > "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_and_total_mapped_read_counts.txt"

# calculate percentages
awk '{print $0, $2/$3}' \
    "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_and_total_mapped_read_counts.txt" \
    > "$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/pver_bams/pver_mapped_read_percents.txt"
```
| Sample                                                                                       | Post-processing_Pver_read_count_condordant_only | Post-processing_total_mapped_read_count | Percent_Pver_reads |
|----------------------------------------------------------------------------------------------|-------------------------------------------------|-----------------------------------------|--------------------|
| 2024_VATI_Pver_21_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 176042818                                       | 336872789                               | 0.52258            |
| 2024_VATI_Pver_22_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 168635188                                       | 328158377                               | 0.513884           |
| 2024_VATI_Pver_23_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 151224676                                       | 292872027                               | 0.516351           |
| 2024_VATI_Pver_24_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 137435262                                       | 271939942                               | 0.505388           |
| 2024_VATI_Pver_25_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 116004972                                       | 229791777                               | 0.504826           |
| 2024_VATI_Pver_26_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 110263344                                       | 216203805                               | 0.509997           |
| 2024_VATI_Pver_27_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 112565296                                       | 222965633                               | 0.504855           |
| 2024_VATI_Pver_28_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 116284158                                       | 235054474                               | 0.494712           |
| 2024_VATI_Pver_29_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 138842314                                       | 272361527                               | 0.509772           |
| 2024_VATI_Pver_30_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 99908544                                        | 195771705                               | 0.510332           |
| 2024_VATI_Pver_31_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 127183126                                       | 252517657                               | 0.50366            |
| 2024_VATI_Pver_32_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 149764782                                       | 299038006                               | 0.500822           |
| 2024_VATI_Pver_33_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 73224156                                        | 144087190                               | 0.508193           |
| 2024_VATI_Pver_34_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 134343792                                       | 268303671                               | 0.500715           |
| 2024_VATI_Pver_35_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 69922392                                        | 137841677                               | 0.507266           |
| 2024_VATI_Pver_36_1_227H3WLT4-L008_bt2_PverCDdedup_primary_minq20_mlen20_pver_reheadered.bam | 72196466                                        | 141689136                               | 0.509541           |


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
#SBATCH --cpus-per-task=30

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

Once genotyping is complete, there will we a separate vcf for each genomic region (scaffold/contig). Combine them into a single multi-sample VCF with BCFtools:
```bash
# make list of vcf files
ls *.vcf.gz | sort > vcf_list.txt


# combine VCFs
crun.bcftools bcftools concat -f vcf_list.txt --threads 34 -Ov -o pver_pilot_combined_genotypes.vcf
```

Now count the total number of called genotypes (variants):
```bash
grep -v "#" pver_pilot_combined_genotypes.vcf | wc -l
```
11,860,400 total SNPs


## 12. Variant filtering

First we need to know a bit about the depth of coverage at each variant
```bash
# Get the total site depth per position, i.e., the sum of read depths from all samples at each variant site, using the INFO/DP field.
crun.bcftools bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' NC_089312.1_Pverrucosa_genotypes.vcf.gz > NC_089312.1_Pverrucosa_genotypes.total_site_depth            # test using one contig
crun.bcftools bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' pver_pilot_combined_genotypes.vcf > pver_pilot_combined_genotypes.total_site_depth                # run on full vcf
```
```
NC_089312.1_Pverrucosa  14681   6
NC_089312.1_Pverrucosa  14695   8
NC_089312.1_Pverrucosa  14719   8
NC_089312.1_Pverrucosa  15847   26
NC_089312.1_Pverrucosa  15864   57
NC_089312.1_Pverrucosa  15869   64
NC_089312.1_Pverrucosa  15889   93
NC_089312.1_Pverrucosa  15901   105
NC_089312.1_Pverrucosa  15960   197
NC_089312.1_Pverrucosa  15966   194
...
```

```bash
# Get the per sample depth at each position
crun.bcftools bcftools query -f '%CHROM\t%POS[\t%DP]\n' NC_089312.1_Pverrucosa_genotypes.vcf.gz > NC_089312.1_Pverrucosa_genotypes.per_sample_depth            # test using one contig
crun.bcftools bcftools query -f '%CHROM\t%POS[\t%DP]\n' pver_pilot_combined_genotypes.vcf > pver_pilot_combined_genotypes.per_sample_depth                # run on full vcf
```
```
NC_089312.1_Pverrucosa  14681   0       0       0       0       2       0       0       1       0       0       3       0       0       0       0       0
NC_089312.1_Pverrucosa  14695   0       0       0       0       2       0       0       1       0       0       3       2       0       0       0       0
NC_089312.1_Pverrucosa  14719   0       0       0       0       0       0       0       3       0       0       3       2       0       0       0       0
```

```
# Get a set of summary stats for the vcf (distributions of quality scores, depths, etc)
crun.bcftools bcftools stats NC_089312.1_Pverrucosa_genotypes.vcf.gz > NC_089312.1_Pverrucosa_genotypes.stats            # test using one contig
crun.bcftools bcftools stats pver_pilot_combined_genotypes.vcf > pver_pilot_combined_genotypes.stats                # run on full vcf
```

Plot depth per site in R:
```r
# SNP depth of coverage analysis - Pver Pilot
# 2025-04-25
# Jason A. Toy

library(tidyverse)

rm(list = ls())
setwd("/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/vcf")


# load total depth file
td <- read_tsv("pver_pilot_combined_genotypes.total_site_depth", col_names = FALSE)

# calculate summary stats
mean_td <- mean(td$X3)
median_td <- median(td$X3)
sd_td <- sd(td$X3)


# plot distribution
p <- ggplot(td) +
  geom_freqpoly(aes(x = X3), binwidth = 1) + 
  xlab("SNP depth") +
  xlim(c(0,5000)) +
  ylab("Count") +
  geom_vline(xintercept = mean_td, color = "blue") +
  geom_vline(xintercept = median_td, color = "darkgreen") +
  geom_vline(xintercept = median_td + sd_td, color = "black", linetype = "dashed") +
  annotate("text", x=3000, y=30000, label=paste0(round(length(td$X3)/10^6,2), " million SNPs before filtering"), color = 'red') +
  annotate("text", x=3000, y=28000, label=paste0("mean total depth = ", round(mean_td, 1)), color = "blue") +
  annotate("text", x=3000, y=26000, label=paste0("median total depth = ", median_td), color = "darkgreen") +
  annotate("text", x=3000, y=24000, label=paste0("--- ", "median total depth + 1 SD = ", round(median_td + sd_td, 1)), color = "black") +
  theme_bw()

ggsave(p, "snp_depth_plot.png", width = "12", height = "8", units = "in")
```
![](snp_depth_plot.png)

### Filter variants with with PLINK2

First filter on quality scores (mapping & variant) and depth with BCFtools:
```bash
crun.bcftools bcftools filter --threads 30 -e 'QUAL < 30 || INFO/MQ < 40 || INFO/DP < 32 || INFO/DP > 1790' pver_pilot_combined_genotypes.vcf -Oz -o pver_pilot_QDPfiltered_genotypes.vcf.gz
```
This leaves **1,603,744** SNPs

Next filter based on missingess and MAF and remove indels and multiallelic SNPs:
```bash
module load plink/2024.03.02    # Note: this is PLINK2

crun.plink plink2 \
  --vcf pver_pilot_QDPfiltered_genotypes.vcf.gz \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno 0.2 \
  --mind 0.2 \
  --maf 0.05 \
  --make-pgen \
  --out pver_pilot_MISSMAFfiltered_genotypes
```

Convert to VCF for viewing:
```bash
crun.plink plink2 \
  --pgen pver_pilot_MISSMAFfiltered_genotypes.pgen \
  --psam pver_pilot_MISSMAFfiltered_genotypes.psam \
  --pvar pver_pilot_MISSMAFfiltered_genotypes.pvar \
  --recode vcf \
  --out pver_pilot_MISSMAFfiltered_genotypes
```
This leaves **1,014,944** SNPs


## 13. LD-pruning with PLINK2
PLINK2 requires unique values in the ID field for each SNP. Our dataset does not contain any ID values, but instead has "." placeholders.
So the first thing we need to do is replace these with a chromosome/position-based ID:
```bash
crun.plink plink2 \
  --pgen pver_pilot_MISSMAFfiltered_genotypes.pgen \
  --psam pver_pilot_MISSMAFfiltered_genotypes.psam \
  --pvar pver_pilot_MISSMAFfiltered_genotypes.pvar \
  --set-all-var-ids @:#:\$r:\$a \
  --make-pgen \
  --out pver_pilot_MISSMAFfiltered_uniqIDs
```
The `@:#:$r:$a` syntax creates IDs formatted as `CHROM:POS:REF:ALT`


Now run LD estimation. Output is an LD-pruned list of SNPs:
```bash
crun.plink plink2 \
  --pgen pver_pilot_MISSMAFfiltered_uniqIDs.pgen \
  --psam pver_pilot_MISSMAFfiltered_uniqIDs.psam \
  --pvar pver_pilot_MISSMAFfiltered_uniqIDs.pvar \
  --indep-pairwise 50 10 0.5 \
  --out pver_pilot_MISSMAFfiltered_ld \
  --bad-ld
```
- 50 → window size in SNPs (can also be specified in kb if 'kb' is added to as a suffix)
- 10 → step size (how many SNPs to shift the window each time)
- 0.2 → r² threshold (SNPs with pairwise r² > 0.2 are considered linked)
- **Note**: I had to add the `--bad-ld` flag for now to force plink2 to run the ld estimation even though sample size is low (<50). This will not be an issue with the full dataset.
- I used r2 > 0.5 to be less stringent to account for errors in LD calculation due to low sample size, but for full dataset, may want to try r2 > 0.2.


Other parameter sets I've seen in the literature:
```
--indep-pairwise 50 10 0.5
--indep-pairwise 50 10 0.1
200kb 20 0.2
--indep-pairwise 50 10 0.2
```


Filter dataset using pruned SNP list:
```bash
crun.plink plink2 \
  --pgen pver_pilot_MISSMAFfiltered_uniqIDs.pgen \
  --pvar pver_pilot_MISSMAFfiltered_uniqIDs.pvar \
  --psam pver_pilot_MISSMAFfiltered_uniqIDs.psam \
  --extract pver_pilot_MISSMAFfiltered_ld.prune.in \
  --make-pgen \
  --out pver_pilot_ld_pruned_genotypes
```

Convert to VCF for viewing/downstream use:
```bash
crun.plink plink2 \
  --pgen pver_pilot_ld_pruned_genotypes.pgen \
  --psam pver_pilot_ld_pruned_genotypes.psam \
  --pvar pver_pilot_ld_pruned_genotypes.pvar \
  --recode vcf \
  --out pver_pilot_ld_pruned_genotypes
```

**85,473** SNPs remaining after LD-pruning


## 13. Principle components analysis

Use PLINK2 to run a PCA:
```bash
crun.plink plink2 \
  --pgen pver_pilot_ld_pruned_genotypes.pgen \
  --psam pver_pilot_ld_pruned_genotypes.psam \
  --pvar pver_pilot_ld_pruned_genotypes.pvar \
  --pca 10 \
  --out pver_pilot_ld_pruned_pca
```
`--pca 10` = calculates the first 10 principal components using exact PCA. For fast approximate PCA (very efficient, but less accurate), use `pca approx 10`.

This gives another error because with a low sample size, PLINK does not trust its allele frequency estimations. So it makes you manually calculate them first if you want to proceed with the PCA anyway:
```bash
crun.plink plink2 \
  --pgen pver_pilot_ld_pruned_genotypes.pgen \
  --psam pver_pilot_ld_pruned_genotypes.psam \
  --pvar pver_pilot_ld_pruned_genotypes.pvar \
  --freq \
  --out pver_pilot_ld_pruned_genotypes
```
This creates a `.afreq` file with calculated allele frequencies.

Now run PCA:
```bash
crun.plink plink2 \
  --pgen pver_pilot_ld_pruned_genotypes.pgen \
  --psam pver_pilot_ld_pruned_genotypes.psam \
  --pvar pver_pilot_ld_pruned_genotypes.pvar \
  --read-freq pver_pilot_ld_pruned_genotypes.afreq \
  --pca 10 \
  --out pver_pilot_ld_pruned_pca
```

Run MDS with PLINK1.9:
```bash
# convert to plink BED format
crun.plink plink2 \
  --pgen pver_pilot_ld_pruned_genotypes.pgen \
  --psam pver_pilot_ld_pruned_genotypes.psam \
  --pvar pver_pilot_ld_pruned_genotypes.pvar \
  --make-bed \
  --out pver_pilot_ld_pruned_genotypes

# switch to PLINK1.9
module load plink/1.9-20240319


# calculate distances and run MDS

# calculate DISSIMILARITY (1-IBS)
crun.plink plink --bfile pver_pilot_ld_pruned_genotypes --distance square 1-ibs flat-missing --out pver_pilot_ld_pruned_genotypes --allow-extra-chr
# calculate similarity (IBS)
crun.plink plink --bfile pver_pilot_ld_pruned_genotypes --distance square ibs flat-missing --out pver_pilot_ld_pruned_genotypes --allow-extra-chr

# --cluster (and therefore --mds-plot) requires a similarity matrix (IBS) not a dissimilarity matrix (1-IBS) so we have to recalculate it:
crun.plink plink --bfile pver_pilot_ld_pruned_genotypes --distance square ibs flat-missing --cluster --mds-plot 10 eigvals --out pver_pilot_ld_pruned_mds --allow-extra-chr
```


Plot PCA and MDS in R:
```r
# Plot PCAs from VCFs - Pver Pilot
# 2025-04-28
# Jason A. Toy

rm(list = ls())

setwd("/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/vcf/")

library(dplyr)
library(ggplot2)
library(ggrepel)

# Load eigenvectors and eigenvalues
eigenvec <- read.delim("pver_pilot_ld_pruned_pca.eigenvec", header = TRUE, sep = "\t")
eigenval <- read.delim("pver_pilot_ld_pruned_pca.eigenval", header = FALSE)

# Calculate percentage of variance explained
percent_var <- (eigenval$V1 / sum(eigenval$V1)) * 100
# 31.910371 19.874622 13.395112 10.226324  9.096963  8.523732  1.946494  1.787711  1.635913  1.602758

# Scree Plot
scree <- percent_var %>% as_tibble %>% rownames_to_column() %>% dplyr::rename(PC = rowname) %>% mutate(PC = as.numeric(PC))

ggplot(scree, aes(x = PC, y = value)) +
  geom_point() + 
  geom_line() + 
  ylab("Percent variance") +
  theme_minimal()
  

# Modify eigenvec for plotting
eigenvec_plot <- eigenvec %>%
  separate(X.IID, into = c("Year", "Location", "Species", "Genotype", "Lib_ID"), sep = "_", extra = "merge") %>% 
  mutate(Geno_ID = paste0(Location, "_", Species, "_", Genotype))


# Plot PCAs

# Just refrence samples (PC 1/2)
PC12 <- ggplot(eigenvec_plot, aes(x = PC1, y = PC2, color = Species)) +
  #scale_color_manual(values = c("#9467BDFF", "#E377C2FF", "#1F77B4FF", "#17BECFFF", "#2CA02CFF", "#BCBD22FF", "#FF7F0EFF", "#D62728FF")) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC1: ", round(percent_var[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2], 2), "% variance")) +
  geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  theme_minimal()

PC34 <- ggplot(eigenvec_plot, aes(x = PC3, y = PC4, color = Species)) +
  #scale_color_manual(values = c("#9467BDFF", "#E377C2FF", "#1F77B4FF", "#17BECFFF", "#2CA02CFF", "#BCBD22FF", "#FF7F0EFF", "#D62728FF")) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC3: ", round(percent_var[3], 2), "% variance")) +
  ylab(paste0("PC4: ", round(percent_var[4], 2), "% variance")) +
  geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  theme_minimal()


# 3D plot, (PCs 1-3)

library(plotly)

plot_ly(
  data = eigenvec_plot,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Species,
  #colors = c("#1F77B4FF", "#2CA02CFF", "#BCBD22FF", "#FF7F0EFF"),
  text = ~Geno_ID,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 3)
) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("PC1: ", round(percent_var[1], 2), "% variance")),
      yaxis = list(title = paste0("PC2: ", round(percent_var[2], 2), "% variance")),
      zaxis = list(title = paste0("PC3: ", round(percent_var[3], 2), "% variance"))
    )
  )


# Plot various PCA combinations side-by-side
library(cowplot)
pdf(file = "pver_pilot_pca_plot.pdf", height = 5, width = 10)
plot_grid(PC12+theme(legend.position = "none"), PC34, ncol = 2)
dev.off()


# Or

library(patchwork)
PC12+theme(legend.position = "none") + PC34



### MDS plot ###

# Load data
mds <- read.delim("pver_pilot_ld_pruned_mds.mds", header = TRUE, sep = "")
mds_eigenval <- read.delim("pver_pilot_ld_pruned_mds.mds.eigvals", header = FALSE)


# Calculate percentage of variance explained
mds_percent_var <- (mds_eigenval$V1 / sum(mds_eigenval$V1)) * 100
# 43.8433102 22.3020291 11.6108199  8.0809120  7.3672850  6.0334479  0.3582680  0.1439877  0.1335114  0.1264289


# Scree Plot
scree <- mds_percent_var %>% as_tibble %>% rownames_to_column() %>% dplyr::rename(Dim = rowname) %>% mutate(Dim = as.numeric(Dim))

ggplot(scree, aes(x = Dim, y = value)) +
  geom_point() + 
  geom_line() + 
  ylab("Percent variance") +
  theme_minimal()


# Modify mds_coords for plotting
mds_plotdata <- mds %>%
  separate(IID, into = c("Year", "Location", "Species", "Genotype", "Lib_ID"), sep = "_", extra = "merge") %>% 
  mutate(Geno_ID = paste0(Location, "_", Species, "_", Genotype))

# Plot MDS
mds12 <- ggplot(mds_plotdata, aes(x = C1, y = C2, color = Species)) +
  #scale_color_manual(values = c("#9467BDFF", "#E377C2FF", "#1F77B4FF", "#17BECFFF", "#2CA02CFF", "#BCBD22FF", "#FF7F0EFF", "#D62728FF")) +
  #ylim(c(-0.125, 0.15)) +
  #xlim(c(-0.2, 0.2)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("MDS1: ", round(mds_percent_var[1], 2), "% variance")) +
  ylab(paste0("MDS2: ", round(mds_percent_var[2], 2), "% variance")) +
  geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = 1000) +
  theme_minimal()

mds34 <- ggplot(mds_plotdata, aes(x = C3, y = C4, color = Species)) +
  #scale_color_manual(values = c("#9467BDFF", "#E377C2FF", "#1F77B4FF", "#17BECFFF", "#2CA02CFF", "#BCBD22FF", "#FF7F0EFF", "#D62728FF")) +
  #ylim(c(-0.125, 0.15)) +
  #xlim(c(-0.2, 0.2)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("MDS3: ", round(mds_percent_var[3], 2), "% variance")) +
  ylab(paste0("MDS4: ", round(mds_percent_var[4], 2), "% variance")) +
  geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = 1000) +
  theme_minimal()


pdf(file = "pver_pilot_mds_plot.pdf", height = 5, width = 10)
plot_grid(mds12+theme(legend.position = "none"), mds34, ncol = 2)
dev.off()

# CONCLUSION: PCA and MDS plots are nearly identical
```
![](pver_pilot_pca_plot.png)
Samples seem to separate clearly into two groups along PC1. We'll see how this pattern holds with the full dataset.

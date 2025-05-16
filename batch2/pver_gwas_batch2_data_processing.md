Processing of 2024 P. verrucosa WGS data from American Samoa using GATK4 tools and adapted best practices - Data Batch 2
================
**AUTHOR:** Jason A. Toy  
**DATE:** 2025-03-27 <br><br>

## Data Prep

Download data, untar each archive, and check md5sums:

```bash
cd /RC/group/rc_barshis_lab/taxonarchive/2024-11-26_TxGenPverrucosa_batch2/

# untar archives
for file in `ls *.tar`; do
  echo $file
  tar -xf $file
  echo "Done extracting $file"
done


cd ./24392Brs_N24193

# check md5sums
md5sum -c MD5SUMS_24392Brs_N24193_L005.txt > md5checks.txt
```


Create sample table (fastq_table) and list (fastq_list):

Sample table columns (in order):

- `prefix` the prefix of raw fastq file names

- `lane_number` lane number; each sequencing lane or batch should be
  assigned a unique identifier. In this case it is "[run_ID]_[lane_id]"

- `seq_id` sequence ID; this variable is relevant when different
  libraries were prepared out of the same sample and were run in the
  same lane (e.g. if you wanted to include a replicate). In this case,
  seq_id should be used to distinguish these separate libraries. If you
  only have a single library prepared from each of your samples (even if
  you sequence that same library across multiple lanes), you can just
  put 1 for all samples in this column. In this case seq_id will refer to
  technical replicates of extracted DNA that we submitted to TxGen. These
  are given values of "1" or "2" (originally "A" and "B" in the submission
  sheet). Replicate libraries prepped by TxGen from the same submitted
  extract are denoted by the `prep_id` column (see below).

- `sample_id` sample ID; a unique identifier for each individual
  sequenced, e.g., "2024_FALU_Pver_01"

- `population` population name; the population or other relevant
  grouping variable that the individual belongs to, e.g., "FALU"

- `data_type` data type; there are only two allowed entries here: pe
  (for paired-end data) or se (for single end data). We need this in the
  table because for some of our processing steps, the commands are
  slightly different for paired-end and single-end data.

-`year` year samples were collected

-`project` project identifier, e.g., ASGWAS

-`location_id` numeric code for sampling location, e.g., "S01"

-`species` species abbreviation of the sampled individual/colony (as identified in the field), e.g., Pver

-`colony_id_by_location` genotype number within a sampling location (01-40)

-`tube_id` genotype number across all sampling locations (tube number)

-`prep_id` identifier for the round of libray preparation during which this library was made. To avoid avoid oversequencing, we asked TxGen to prep two separate libraries for each sample ("A" and "B"). Some samples failed in one prep but not the other and were reprepped in a third round, "C".

-`seq_project_id` vendor (TxGen)-assigned project number

-`run_id` ID assigned to the sequencing run by TxGen

-`lane_id` lane ID within a sequencing run (in this case L001 - L008)

-`vendor_id` ID assigned to the library by TxGen, e.g., "S421"



First few lines of fastq_table_pver_gwas_batch2.txt:
```
prefix  lane_number     seq_id  sample_id       population      data_type       year    project location_id     species colony_id_by_location   tube_id prep_id seq_project_id  run_id  lane_id vendor_id
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24193_S421_L007 R24193_L007     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     B       24392Brs        R24193  L007    S421
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24193_S421_L008 R24193_L008     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     B       24392Brs        R24193  L008    S421
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24193_S25_L005    R24193_L005     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     A       24392Brs        R24193  L005    S25
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24193_S25_L006    R24193_L006     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     A       24392Brs        R24193  L006    S25
24392Brs_2024-ASGWAS-S01-Pver-02-542-B_R24193_S422_L007 R24193_L007     1       2024_FALU_Pver_02       FALU    pe      2024    ASGWAS  S01     Pver    02      542     B       24392Brs        R24193  L007    S422
24392Brs_2024-ASGWAS-S01-Pver-02-542-B_R24193_S422_L008 R24193_L008     1       2024_FALU_Pver_02       FALU    pe      2024    ASGWAS  S01     Pver    02      542     B       24392Brs        R24193  L008    S422
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24193_S26_L005    R24193_L005     1       2024_FALU_Pver_02       FALU    pe      2024    ASGWAS  S01     Pver    02      542     A       24392Brs        R24193  L005    S26
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24193_S26_L006    R24193_L006     1       2024_FALU_Pver_02       FALU    pe      2024    ASGWAS  S01     Pver    02      542     A       24392Brs        R24193  L006    S26
24392Brs_2024-ASGWAS-S01-Pver-03-543-B_R24193_S423_L007 R24193_L007     1       2024_FALU_Pver_03       FALU    pe      2024    ASGWAS  S01     Pver    03      543     B       24392Brs        R24193  L007    S423
```

Make fastq_list:
``` bash
cut -f1 fastq_table_pver_gwas_batch2.txt | tail -n +2 > fastq_list_pver_gwas_batch2.txt
```



Assign starting environment variables and check environment and fastq_table:

``` bash
BASEDIR=/archive/barshis/barshislab/jtoy/pver_gwas/
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_pver_gwas_batch2.txt
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table_pver_gwas_batch2.txt



#Use a for loop to print the sample IDs for all sample files in your sample list:

for SAMPLEFILE in `cat $SAMPLELIST`; do   # Loop through each of the prefixes listed in our fastq list

    # For each prefix, extract the associated sample ID (column 4) and population (column 5) from the table
    SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
    POPULATION=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
    echo $SAMPLEFILE refers to sample $SAMPLE_ID from $POPULATION

done
```

```
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24193_S421_L007 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24193_S421_L008 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24193_S25_L005 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24193_S25_L006 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-02-542-B_R24193_S422_L007 refers to sample 2024_FALU_Pver_02 from FALU
24392Brs_2024-ASGWAS-S01-Pver-02-542-B_R24193_S422_L008 refers to sample 2024_FALU_Pver_02 from FALU
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24193_S26_L005 refers to sample 2024_FALU_Pver_02 from FALU
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24193_S26_L006 refers to sample 2024_FALU_Pver_02 from FALU
24392Brs_2024-ASGWAS-S01-Pver-03-543-B_R24193_S423_L007 refers to sample 2024_FALU_Pver_03 from FALU
24392Brs_2024-ASGWAS-S01-Pver-03-543-B_R24193_S423_L008 refers to sample 2024_FALU_Pver_03 from FALU
...
```



## Check read quality and trim reads using fastp (v0.23.2)

``` bash
mkdir $BASEDIR/reports/fastp_reports
```

 

Create array SLURM script to run fastp on all samples at once:

``` bash
nano fastp_array.slurm
```

``` bash
#!/bin/bash

#SBATCH --job-name fastp_array_2024-03-27
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-1560%16
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 7-00:00:00


## Load modules
module load container_env
module load fastp

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
RAWDATA=/RC/group/rc_barshis_lab/taxonarchive/2024-11-26_TxGenPverrucosa_batch2/24392Brs_N24193 #path to raw fq.gz files
OUTDIR=$BASEDIR/pver_gwas/trimmed_fastq
SAMPLELIST=$BASEDIR/pver_gwas/sample_lists/fastq_list_pver_gwas_batch2.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table (without the header line).
SAMPLETABLE=$BASEDIR/pver_gwas/sample_lists/fastq_table_pver_gwas_batch2.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
RAW_R1=_R1_001.fastq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.
RAW_R2=_R2_001.fastq.gz # Suffix to raw fastq files. Use reverse reads with paired-end data.


## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library. This is for the naming of trimmed/processed files
SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
PREP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 13`
SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$PREP_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely

echo $SAMPLE_UNIQ_ID

## Define the output path and file prefix
SAMPLETRIM=$OUTDIR/$SAMPLE_UNIQ_ID


## Run fastp
crun.fastp fastp -i $RAWDATA/$SAMPLEFILE$RAW_R1 \
    -I $RAWDATA/$SAMPLEFILE$RAW_R2 \
    -o ${SAMPLETRIM}_f_paired_trim.fastq.gz \
    -O ${SAMPLETRIM}_r_paired_trim.fastq.gz \
    --adapter_fasta $BASEDIR/pver_gwas/adapters.fa \
    --cut_tail \
    --trim_poly_g \
    -l 40 \
    -h $BASEDIR/pver_gwas/reports/fastp_reports/${SAMPLE_UNIQ_ID}_fastp.html \
    -j $BASEDIR/pver_gwas/reports/fastp_reports/${SAMPLE_UNIQ_ID}_fastp.json \
    --thread 16

```

``` bash
sbatch fastp_array.slurm
```

 

Move html reports to new subdirectory:

``` bash
mkdir $BASEDIR/reports/fastp_reports/html_reports

mv *.html html_reports
```

 

### Visualize with MultiQC (v1.13)

Run multiqc on the directory with the json report files:

``` bash
module load container_env
module load multiqc/1.13

crun.multiqc multiqc --interactive --filename multiqc_report_pver_batch2_fastp_interactive .
# --interactive forces the creation of an interactive html file even when sample sizes are high (instead of a flat file)
```


The following files had fewer than 100,000 reads:

| File                                                             | Passed Filter | Low Quality | Too Many N | Too Short | Too Long | Percent Passed Filter |
|------------------------------------------------------------------|---------------|-------------|------------|-----------|----------|-----------------------|
| 24392Brs_2024-ASGWAS-S10-Pver-03-823_R24193_S310_L005_R1_001     | 0             | 52          | 0          | 256       | 0        | 0.00                  |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823_R24193_S310_L006_R1_001     | 2             | 52          | 0          | 266       | 0        | 0.63                  |
| 24392Brs_2024-ASGWAS-S06-Pver-14-614_R24193_S99_L006_R1_001      | 6             | 0           | 0          | 0         | 0        | 100.00                |
| 24392Brs_2024-ASGWAS-S02-Pver-11-871-B_R24193_S754_L008_R1_001   | 10            | 0           | 0          | 0         | 0        | 100.00                |
| 24392Brs_2024-ASGWAS-S02-Pver-11-871-B_R24193_S754_L007_R1_001   | 12            | 0           | 0          | 0         | 0        | 100.00                |
| 24392Brs_2024-ASGWAS-S06-Pver-14-614_R24193_S99_L005_R1_001      | 20            | 0           | 0          | 2         | 0        | 90.91                 |
| 24392Brs_2024-ASGWAS-S06-Pver-24-624_R24193_S109_L006_R1_001     | 184           | 2           | 0          | 0         | 0        | 98.92                 |
| 24392Brs_2024-ASGWAS-S06-Pver-24-624_R24193_S109_L005_R1_001     | 194           | 2           | 0          | 0         | 0        | 98.98                 |
| 24392Brs_2024-ASGWAS-S06-Pver-27-627_R24193_S112_L006_R1_001     | 436           | 0           | 2          | 0         | 0        | 99.54                 |
| 24392Brs_2024-ASGWAS-S06-Pver-27-627_R24193_S112_L005_R1_001     | 438           | 8           | 0          | 4         | 0        | 97.33                 |
| 24392Brs_2024-ASGWAS-S06-Pver-16-616_R24193_S101_L006_R1_001     | 520           | 2           | 0          | 10        | 0        | 97.74                 |
| 24392Brs_2024-ASGWAS-S06-Pver-16-616_R24193_S101_L005_R1_001     | 564           | 18          | 0          | 12        | 0        | 94.95                 |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823-B_R24193_S706_L007_R1_001   | 726           | 66538       | 0          | 384924    | 0        | 0.16                  |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823-B_R24193_S706_L008_R1_001   | 730           | 67944       | 4          | 370340    | 0        | 0.17                  |
| 24392Brs_2024-ASGWAS-S02-Ahya-21-915_R24193_S369_L006_R1_001     | 2468          | 16          | 0          | 4         | 0        | 99.20                 |
| 24392Brs_2024-ASGWAS-S02-Ahya-21-915_R24193_S369_L005_R1_001     | 2528          | 24          | 0          | 8         | 0        | 98.75                 |
| 24392Brs_2024-ASGWAS-S11-Pspp-Extra1-918_R24193_S372_L006_R1_001 | 3952          | 16          | 0          | 14        | 0        | 99.25                 |
| 24392Brs_2024-ASGWAS-S11-Pspp-Extra1-918_R24193_S372_L005_R1_001 | 4092          | 36          | 0          | 20        | 0        | 98.65                 |
| 24392Brs_2024-ASGWAS-S06-Pver-12-612_R24193_S97_L005_R1_001      | 9586          | 272         | 0          | 546       | 0        | 92.14                 |
| 24392Brs_2024-ASGWAS-S06-Pver-12-612_R24193_S97_L006_R1_001      | 9766          | 282         | 0          | 442       | 0        | 93.10                 |
| 24392Brs_2024-ASGWAS-S06-Pver-18-618_R24193_S103_L005_R1_001     | 34056         | 802         | 4          | 556       | 0        | 96.15                 |
| 24392Brs_2024-ASGWAS-S06-Pver-18-618_R24193_S103_L006_R1_001     | 34800         | 852         | 0          | 534       | 0        | 96.17                 |


Sample 823 (S10_Pver_03) was clearly problematic.
But overall, sequences look good:
![image](https://github.com/user-attachments/assets/14f6689d-4d7a-47a8-9af7-52d2780b39ac)

![image](https://github.com/user-attachments/assets/e8753fc5-dd65-4391-9c3a-93fffa119891)

<br>

## Map reads to hologenome

`hologenome_mapping_array.slurm`
```bash
#!/bin/bash

#SBATCH --job-name hologenome_mapping_array_pverbatch2_2025-05-09
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-1560%100
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 7-00:00:00


## Load modules
module load container_env bowtie2

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
FASTQDIR=$BASEDIR/pver_gwas/pver_gwas_batch2/trimmed_fastq #path to trimmed fastq.gz files
OUTDIR=$BASEDIR/pver_gwas/pver_gwas_batch2/bam
SAMPLELIST=$BASEDIR/pver_gwas/pver_gwas_batch2/sample_lists/fastq_list_pver_gwas_batch2.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/pver_gwas/pver_gwas_batch2/sample_lists/fastq_table_pver_gwas_batch2.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns has to be unique. The 6th column should be data type, which is either pe or se.
FASTQ_SUFFIX_1=_f_paired_trim.fastq.gz # Suffix to trimmed fastq files. Use forward reads with paired-end data.
FASTQ_SUFFIX_2=_r_paired_trim.fastq.gz # Suffix to trimmed fastq files. Use reverse reads with paired-end data.
REFBASENAME=combined_pver_cd_hologenome

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library. This is for the naming of trimmed/processed files
SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
PREP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 13`
SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$PREP_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2` # Define platform unit (PU), which is the lane number

echo $SAMPLE_UNIQ_ID

## Define the output path and file prefix
SAMPLEOUT=$OUTDIR/$SAMPLE_UNIQ_ID


## Run bowtie2 to map paired-end reads
# Note: we ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.
crun.bowtie2 bowtie2 -q --phred33 --very-sensitive -p 16 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA \
        -x /cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/$REFBASENAME -1 $FASTQDIR/$SAMPLE_UNIQ_ID$FASTQ_SUFFIX_1 -2 $FASTQDIR/$SAMPLE_UNIQ_ID$FASTQ_SUFFIX_2 \
        -S $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'

#add version number to bowtie2 command header line in SAM file
sed 's/PN:bowtie2\tVN:/PN:bowtie2\tVN:2.4.1/' $SAMPLEOUT'_bt2_'$REFBASENAME'.sam' > $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered.sam'

#remove SAM file
rm $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'

## Change modules
module unload bowtie2
module load container_env gatk
GATK='crun.gatk gatk'

## Query-sort for duplicate removal with GATK
# Run SortSam to sort by query name and convert to BAM
$GATK --java-options "-Xmx100G" SortSam \
  --INPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered.sam' \
  --OUTPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted.bam' \
  --SORT_ORDER queryname

# Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted.val' \
  -M VERBOSE

#remove SAM file
rm $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered.sam'

## Mark and remove duplicates
$GATK --java-options "-Xmx100G" MarkDuplicates \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dedup.bam' \
  --METRICS_FILE $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dupstat.txt' \
  --REMOVE_DUPLICATES true

## Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dedup.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dedup.val' \
  -M VERBOSE

#remove old BAM file
rm $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted.bam'

## Run SortSam to sort by coordinate for downstream processing
$GATK --java-options "-Xmx100G" SortSam \
  --INPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dedup.bam' \
  --OUTPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dedup_coordsorted.bam' \
  --SORT_ORDER coordinate

# Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dedup_coordsorted.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dedup_coordsorted.val' \
  -M VERBOSE

#remove old BAM files
rm $SAMPLEOUT'_bt2_'$REFBASENAME'_reheadered_qsorted_dedup.bam'
```


Run the `parse_bowtie2_output.sh` script to parse summary mapping statistics from each file. The output file is called bowtie_mapping_summary.tsv
```bash
#!/bin/bash

ls *hologenome_mapping*.err > errfile_list.txt

# create data file and add header if it is empty
header="jobid\tsample\tconcordantly_0_times\tconcordantly_1_time\tconcordantly_2_or_more_times\toverall_rate" #assign header value
outfile="bowtie_mapping_summary.tsv"

if [ ! -s "$outfile" ]; then
  # file is empty or does not exist
    echo -e "$header" > "$outfile"
fi

for FILE in `cat errfile_list.txt`; do
    # parse mapping rates from botwtie2 otuput
    jobid=$(echo $FILE | cut -d "_" -f1,2)
    sample=$(head -3 $jobid'_hologenome_mapping_array_pverbatch2_2025-05-09.out' | tail -1)
    con0=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pverbatch2_2025-05-09.err' | head -2 | tail -1)
    con1=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pverbatch2_2025-05-09.err' | head -3 | tail -1)
    con2=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pverbatch2_2025-05-09.err' | head -4 | tail -1)
    overall=$(grep -oP '\d+\.\d+%? overall alignment rate' $jobid'_hologenome_mapping_array_pverbatch2_2025-05-09.err' | cut -d" " -f1)

    # Append data to output file
    echo -e "$jobid\t$sample\t$con0\t$con1\t$con2\t$overall" >> "$outfile"

done
```

Check for low-mapping samples:
```bash
sort -V -k6,6 bowtie_mapping_summary.tsv | grep -v "Ahya" | less -S
```

```
jobid   sample  concordantly_0_times    concordantly_1_time     concordantly_2_or_more_times    overall_rate
4419996_1403    2024_AOAA_Pver_03_1_A_R24193_L005       4.5     0.0     4.5     0.00%
4419996_1404    2024_AOAA_Pver_03_1_A_R24193_L006       100.00% 0.00%   0.00%   0.00%
4419996_1402    2024_AOAA_Pver_03_1_B_R24193_L008       100.00% 0.00%   0.00%   4.25%
4419996_1401    2024_AOAA_Pver_03_1_B_R24193_L007       100.00% 0.00%   0.00%   6.06%
4419996_747     2024_OFU3_Pver_13_1_B_R24193_L007       47.16%  34.22%  18.61%  63.48%
4419996_748     2024_OFU3_Pver_13_1_B_R24193_L008       47.24%  34.17%  18.59%  63.49%
4419996_749     2024_OFU3_Pver_13_1_A_R24193_L005       44.75%  35.60%  19.65%  64.20%
4419996_750     2024_OFU3_Pver_13_1_A_R24193_L006       44.97%  35.53%  19.50%  64.20%
4419996_579     2024_MALO_Pver_11_1_B_R24193_L007       38.95%  39.25%  21.80%  69.49%
4419996_580     2024_MALO_Pver_11_1_B_R24193_L008       38.97%  39.25%  21.78%  69.51%
4419996_581     2024_MALO_Pver_11_1_A_R24193_L005       37.92%  39.77%  22.31%  70.24%
4419996_582     2024_MALO_Pver_11_1_A_R24193_L006       38.16%  39.67%  22.17%  70.24%
4419996_1177    2024_FTEL_Pver_33_1_B_R24193_L007       40.38%  39.13%  20.49%  73.68%
4419996_1178    2024_FTEL_Pver_33_1_B_R24193_L008       40.48%  39.06%  20.47%  73.68%
4419996_1180    2024_FTEL_Pver_33_1_A_R24193_L006       39.18%  39.55%  21.27%  73.76%
4419996_1179    2024_FTEL_Pver_33_1_A_R24193_L005       38.85%  39.69%  21.46%  73.78%
4419996_1548    2024_AOAA_Pver_39_1_B_R24193_L008       37.75%  41.42%  20.83%  76.02%
4419996_1547    2024_AOAA_Pver_39_1_B_R24193_L007       37.62%  41.50%  20.88%  76.05%
4419996_1550    2024_AOAA_Pver_39_1_A_R24193_L006       36.00%  42.35%  21.65%  76.40%
```

Most libraries had 70-90% mapping rate. The only libraries with really low mapping were the problematic AOAA_Pver_03 libraries (sample 823), which had little to no sequence data.

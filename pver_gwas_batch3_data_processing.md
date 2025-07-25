Processing of 2024 P. verrucosa WGS data from American Samoa using GATK4 tools and adapted best practices - Data Batch 3
================
**AUTHOR:** Jason A. Toy  
**DATE:** 2025-05-05 <br><br>

## Data Prep

Download data, untar each archive, and check md5sums:

```bash
cd /RC/group/rc_barshis_lab/taxonarchive/2024-12-05_TxGenPverrucosa_batch3/

# untar archives
for file in `ls *.tar`; do
  echo $file
  tar -xf $file
  echo "Done extracting $file"
done


cd ./24392Brs_N24196

# check md5sums
md5sum -c MD5SUMS_24392Brs_N24196_L001.txt > md5checks.txt

for FILE in MD5SUMS_24392Brs_N24196_L*.txt; do
  md5sum -c $FILE >> md5checks.txt
done
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



First few lines of fastq_table_pver_gwas_batch3.txt:
```
prefix  lane_number     seq_id  sample_id       population      data_type       year    project location_id     species colony_id_by_location   tube_id prep_id seq_project_id  run_id  lane_id vendor_id
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24196_S25_L001    R24196_L001     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     A       24392Brs        R24196  L001    S25
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24196_S25_L002    R24196_L002     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     A       24392Brs        R24196  L002    S25
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24196_S25_L003    R24196_L003     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     A       24392Brs        R24196  L003    S25
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24196_S25_L004    R24196_L004     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     A       24392Brs        R24196  L004    S25
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24196_S421_L005 R24196_L005     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     B       24392Brs        R24196  L005    S421
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24196_S421_L006 R24196_L006     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     B       24392Brs        R24196  L006    S421
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24196_S421_L007 R24196_L007     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     B       24392Brs        R24196  L007    S421
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24196_S421_L008 R24196_L008     1       2024_FALU_Pver_01       FALU    pe      2024    ASGWAS  S01     Pver    01      541     B       24392Brs        R24196  L008    S421
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24196_S26_L001    R24196_L001     1       2024_FALU_Pver_02       FALU    pe      2024    ASGWAS  S01     Pver    02      542     A       24392Brs        R24196  L001    S26
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24196_S26_L002    R24196_L002     1       2024_FALU_Pver_02       FALU    pe      2024    ASGWAS  S01     Pver    02      542     A       24392Brs        R24196  L002    S26
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24196_S26_L003    R24196_L003     1       2024_FALU_Pver_02       FALU    pe      2024    ASGWAS  S01     Pver    02      542     A       24392Brs        R24196  L003    S26
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24196_S26_L004    R24196_L004     1       2024_FALU_Pver_02       FALU    pe      2024    ASGWAS  S01     Pver    02      542     A       24392Brs        R24196  L004    S26
```

Make fastq_list:
``` bash
cut -f1 fastq_table_pver_gwas_batch3.txt | tail -n +2 > fastq_list_pver_gwas_batch3.txt
```


Assign starting environment variables and check environment and fastq_table:

``` bash
BASEDIR=/archive/barshis/barshislab/jtoy/pver_gwas/pver_gwas_batch3
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_pver_gwas_batch3.txt
SAMPLETABLE=$BASEDIR/sample_lists/fastq_table_pver_gwas_batch3.txt



#Use a for loop to print the sample IDs for all sample files in your sample list:

for SAMPLEFILE in `cat $SAMPLELIST`; do   # Loop through each of the prefixes listed in our fastq list

    # For each prefix, extract the associated sample ID (column 4) and population (column 5) from the table
    SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
    POPULATION=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
    echo $SAMPLEFILE refers to sample $SAMPLE_ID from $POPULATION

done
```

```
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24196_S25_L001 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24196_S25_L002 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24196_S25_L003 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541_R24196_S25_L004 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24196_S421_L005 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24196_S421_L006 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24196_S421_L007 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-01-541-B_R24196_S421_L008 refers to sample 2024_FALU_Pver_01 from FALU
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24196_S26_L001 refers to sample 2024_FALU_Pver_02 from FALU
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24196_S26_L002 refers to sample 2024_FALU_Pver_02 from FALU
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24196_S26_L003 refers to sample 2024_FALU_Pver_02 from FALU
24392Brs_2024-ASGWAS-S01-Pver-02-542_R24196_S26_L004 refers to sample 2024_FALU_Pver_02 from FALU
...
```


Check number of sequence file pairs:
```bash
tail -n +2 fastq_table_pver_gwas_batch3.txt | wc -l
```
```
3118
```

Check number of unique libraries:
```bash
cut -f4,3,13 fastq_table_pver_gwas_batch3.txt | tail -n +2 | uniq | wc -l
```
```
780
```

Count number of "A" libraries:
```bash
# counting separate extractions
grep -P "\tA\t" fastq_table_pver_gwas_batch3.txt | cut -f1 | cut -f2 -d "_" | sort | uniq | wc -l

# merging replicate extractions (e.g., 2024-ASGWAS-S02-Pver-14-874-a & 2024-ASGWAS-S02-Pver-14-874-b)
grep -P "\tA\t" fastq_table_pver_gwas_batch3.txt | cut -f1 | cut -f2 -d "_" | cut -f1-6 -d "-" | sort | uniq | wc -l
```
```
384  # counting each extraction separately
380  # excluding replicate extractions
```

Count number of "B" libraries:
```bash
# counting separate extractions
grep -P "\tB\t" fastq_table_pver_gwas_batch3.txt | cut -f1 | cut -f2 -d "_" | sort | uniq | wc -l

# merging replicate extractions
grep -P "\tB\t" fastq_table_pver_gwas_batch3.txt | cut -f1 | cut -f2 -d "_" | sort | uniq | grep -v -P "\-b\-" | wc -l
```
```
382  # counting each extraction separately
378  # excluding replicate extractions
```

Count number of "C" libraries:
```
# counting separate extractions
grep -P "\tC\t" fastq_table_pver_gwas_batch3.txt | cut -f1 | cut -f2 -d "_" | sort | uniq | wc -l

# merging replicate extractions
grep -P "\tC\t" fastq_table_pver_gwas_batch3.txt | cut -f1 | cut -f2 -d "_" | sort | uniq | grep -v -P "\-b\-" | wc -l
```
```
14  # counting each extraction separately
14  # excluding replicate extractions
```


Check number of unique extractions:
```bash
cut -f4,3 fastq_table_pver_gwas_batch3.txt | tail -n +2 | uniq | wc -l
```
```
384
```
We submitted 4 plates of AS samples (including a few some Ahya from Leone), so this matches expectations.


Check number of unique samples:
```bash
cut -f4 fastq_table_pver_gwas_batch3.txt | tail -n +2 | sort | uniq | wc -l
```
```
380
```
We submitted 4 samples with two replicate extractions, so this matches expectations.





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

#SBATCH --job-name fastp_array_2025-05-06
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-3118%100
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 7-00:00:00


## Load modules
module load container_env
module load fastp

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
RAWDATA=/RC/group/rc_barshis_lab/taxonarchive/2024-12-05_TxGenPverrucosa_batch3/24392Brs_N24196 #path to raw fq.gz files
OUTDIR=$BASEDIR/pver_gwas/pver_gwas_batch3/trimmed_fastq
SAMPLELIST=$BASEDIR/pver_gwas/pver_gwas_batch3/sample_lists/fastq_list_pver_gwas_batch3.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table (without the header line).
SAMPLETABLE=$BASEDIR/pver_gwas/pver_gwas_batch3/sample_lists/fastq_table_pver_gwas_batch3.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
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
    --adapter_fasta $BASEDIR/pver_gwas/pver_gwas_batch3/adapters.fa \
    --cut_tail \
    --trim_poly_g \
    -l 40 \
    -h $BASEDIR/pver_gwas/pver_gwas_batch3/reports/fastp_reports/${SAMPLE_UNIQ_ID}_fastp.html \
    -j $BASEDIR/pver_gwas/pver_gwas_batch3/reports/fastp_reports/${SAMPLE_UNIQ_ID}_fastp.json \
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

crun.multiqc multiqc --interactive --filename multiqc_report_pver_batch3_fastp_interactive .
# --interactive forces the creation of an interactive html file even when sample sizes are high (instead of a flat file)
```

The following files had fewer than 100,000 reads:
<br>
| File                                                             | Passed Filter | Low Quality | Too Many N | Too Short | Too Long | Percent Passed Filter |
|------------------------------------------------------------------|---------------|-------------|------------|-----------|----------|-----------------------|
| 24392Brs_2024-ASGWAS-S10-Pver-03-823_R24196_S310_L003_R1_001     | 0             | 78          | 0          | 288       | 0        | 0                     |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823_R24196_S310_L004_R1_001     | 0             | 58          | 0          | 296       | 0        | 0                     |
| 24392Brs_2024-ASGWAS-S02-Pver-11-871-B_R24196_S754_L005_R1_001   | 2             | 0           | 0          | 0         | 0        | 100                   |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823_R24196_S310_L001_R1_001     | 2             | 66          | 0          | 310       | 0        | 0.529101              |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823_R24196_S310_L002_R1_001     | 2             | 50          | 0          | 332       | 0        | 0.520833              |
| 24392Brs_2024-ASGWAS-S06-Pver-14-614_R24196_S99_L002_R1_001      | 4             | 0           | 0          | 0         | 0        | 100                   |
| 24392Brs_2024-ASGWAS-S06-Pver-14-614_R24196_S99_L001_R1_001      | 6             | 0           | 0          | 0         | 0        | 100                   |
| 24392Brs_2024-ASGWAS-S06-Pver-14-614_R24196_S99_L004_R1_001      | 12            | 0           | 0          | 0         | 0        | 100                   |
| 24392Brs_2024-ASGWAS-S06-Pver-14-614_R24196_S99_L003_R1_001      | 14            | 0           | 0          | 0         | 0        | 100                   |
| 24392Brs_2024-ASGWAS-S02-Pver-11-871-B_R24196_S754_L006_R1_001   | 22            | 0           | 0          | 0         | 0        | 100                   |
| 24392Brs_2024-ASGWAS-S06-Pver-24-624_R24196_S109_L003_R1_001     | 138           | 2           | 0          | 0         | 0        | 98.57143              |
| 24392Brs_2024-ASGWAS-S06-Pver-24-624_R24196_S109_L002_R1_001     | 142           | 6           | 0          | 0         | 0        | 95.94595              |
| 24392Brs_2024-ASGWAS-S06-Pver-24-624_R24196_S109_L004_R1_001     | 162           | 0           | 0          | 4         | 0        | 97.59036              |
| 24392Brs_2024-ASGWAS-S06-Pver-24-624_R24196_S109_L001_R1_001     | 184           | 6           | 0          | 0         | 0        | 96.84211              |
| 24392Brs_2024-ASGWAS-S06-Pver-27-627_R24196_S112_L003_R1_001     | 352           | 8           | 0          | 10        | 0        | 95.13514              |
| 24392Brs_2024-ASGWAS-S06-Pver-27-627_R24196_S112_L001_R1_001     | 422           | 14          | 0          | 10        | 0        | 94.61883              |
| 24392Brs_2024-ASGWAS-S06-Pver-27-627_R24196_S112_L004_R1_001     | 472           | 4           | 0          | 6         | 0        | 97.92531              |
| 24392Brs_2024-ASGWAS-S06-Pver-16-616_R24196_S101_L002_R1_001     | 504           | 24          | 0          | 18        | 0        | 92.30769              |
| 24392Brs_2024-ASGWAS-S06-Pver-16-616_R24196_S101_L004_R1_001     | 508           | 2           | 0          | 18        | 0        | 96.21212              |
| 24392Brs_2024-ASGWAS-S06-Pver-27-627_R24196_S112_L002_R1_001     | 522           | 10          | 0          | 4         | 0        | 97.38806              |
| 24392Brs_2024-ASGWAS-S06-Pver-16-616_R24196_S101_L003_R1_001     | 570           | 14          | 0          | 18        | 0        | 94.68439              |
| 24392Brs_2024-ASGWAS-S06-Pver-16-616_R24196_S101_L001_R1_001     | 608           | 10          | 0          | 6         | 0        | 97.4359               |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823-B_R24196_S706_L005_R1_001   | 694           | 46464       | 0          | 399394    | 0        | 0.155413              |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823-B_R24196_S706_L007_R1_001   | 700           | 46918       | 0          | 406316    | 0        | 0.154207              |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823-B_R24196_S706_L006_R1_001   | 702           | 47120       | 0          | 403108    | 0        | 0.155678              |
| 24392Brs_2024-ASGWAS-S10-Pver-03-823-B_R24196_S706_L008_R1_001   | 710           | 52900       | 0          | 405188    | 0        | 0.154752              |
| 24392Brs_2024-ASGWAS-S02-Ahya-21-915_R24196_S369_L003_R1_001     | 2564          | 24          | 0          | 14        | 0        | 98.53958              |
| 24392Brs_2024-ASGWAS-S02-Ahya-21-915_R24196_S369_L002_R1_001     | 2626          | 30          | 0          | 4         | 0        | 98.7218               |
| 24392Brs_2024-ASGWAS-S02-Ahya-21-915_R24196_S369_L004_R1_001     | 2632          | 30          | 0          | 10        | 0        | 98.50299              |
| 24392Brs_2024-ASGWAS-S02-Ahya-21-915_R24196_S369_L001_R1_001     | 2758          | 16          | 4          | 6         | 0        | 99.06609              |
| 24392Brs_2024-ASGWAS-S11-Pspp-Extra1-918_R24196_S372_L003_R1_001 | 4044          | 36          | 0          | 8         | 0        | 98.92368              |
| 24392Brs_2024-ASGWAS-S11-Pspp-Extra1-918_R24196_S372_L004_R1_001 | 4128          | 44          | 0          | 14        | 0        | 98.61443              |
| 24392Brs_2024-ASGWAS-S11-Pspp-Extra1-918_R24196_S372_L002_R1_001 | 4138          | 28          | 0          | 14        | 0        | 98.99522              |
| 24392Brs_2024-ASGWAS-S11-Pspp-Extra1-918_R24196_S372_L001_R1_001 | 4170          | 14          | 0          | 8         | 0        | 99.47519              |
| 24392Brs_2024-ASGWAS-S06-Pver-12-612_R24196_S97_L004_R1_001      | 9612          | 282         | 0          | 536       | 0        | 92.15724              |
| 24392Brs_2024-ASGWAS-S06-Pver-12-612_R24196_S97_L002_R1_001      | 9636          | 250         | 0          | 542       | 0        | 92.40506              |
| 24392Brs_2024-ASGWAS-S06-Pver-12-612_R24196_S97_L003_R1_001      | 9668          | 272         | 0          | 514       | 0        | 92.48135              |
| 24392Brs_2024-ASGWAS-S06-Pver-12-612_R24196_S97_L001_R1_001      | 9758          | 282         | 4          | 558       | 0        | 92.03924              |
| 24392Brs_2024-ASGWAS-S06-Pver-18-618_R24196_S103_L003_R1_001     | 34428         | 742         | 16         | 552       | 0        | 96.33443              |
| 24392Brs_2024-ASGWAS-S06-Pver-18-618_R24196_S103_L002_R1_001     | 34884         | 866         | 4          | 602       | 0        | 95.95115              |
| 24392Brs_2024-ASGWAS-S06-Pver-18-618_R24196_S103_L004_R1_001     | 34898         | 818         | 8          | 540       | 0        | 96.23318              |
| 24392Brs_2024-ASGWAS-S06-Pver-18-618_R24196_S103_L001_R1_001     | 36444         | 822         | 2          | 588       | 0        | 96.27008              |

<br>

Sample 823 (S10_Pver_03) was again problematic, but again, most samples looked good:

<br>

![image](https://github.com/user-attachments/assets/71d471f6-3f61-4fb7-a5b1-b1bba96d864a)
![image](https://github.com/user-attachments/assets/d5e4279b-9138-4fd2-ac57-ed5fcfb04736)

<br>

## Map reads to hologenome

`hologenome_mapping_array.slurm`
```bash
#!/bin/bash

#SBATCH --job-name hologenome_mapping_array_pverbatch3_2025-05-09
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-3118%110
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 14-00:00:00


## Load modules
module load container_env bowtie2

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
FASTQDIR=$BASEDIR/pver_gwas/pver_gwas_batch3/trimmed_fastq #path to trimmed fastq.gz files
OUTDIR=$BASEDIR/pver_gwas/pver_gwas_batch3/bam
SAMPLELIST=$BASEDIR/pver_gwas/pver_gwas_batch3/sample_lists/fastq_list_pver_gwas_batch3.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/pver_gwas/pver_gwas_batch3/sample_lists/fastq_table_pver_gwas_batch3.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns has to be unique. The 6th column should be data type, which is either pe or se.
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
This took about **10 days** to complete running 110 array jobs simultaneously. Only one job failed: `4432526_2104`, corresponding to `24392Brs_2024-ASGWAS-S08-Pver-07-707-B_R24196_S589_L006`. It appears that it ran out of memory during the deduplication step, even though 100GB was allocated. I reran this sample alone, this time allocating 120GB and using 38 threads for the bowtie step instead of 16. This script is called `hologenome_mapping_sample707redo.slurm`.

<br>

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
    sample=$(head -3 $jobid'_hologenome_mapping_array_pverbatch3_2025-05-15.out' | tail -1)
    con0=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pverbatch3_2025-05-15.err' | head -2 | tail -1)
    con1=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pverbatch3_2025-05-15.err' | head -3 | tail -1)
    con2=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pverbatch3_2025-05-15.err' | head -4 | tail -1)
    overall=$(grep -oP '\d+\.\d+%? overall alignment rate' $jobid'_hologenome_mapping_array_pverbatch3_2025-05-15.err' | cut -d" " -f1)

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
4432526_2775    2024_AOAA_Pver_03_1_A_R24196_L001       100.00% 0.00%   0.00%   0.00%
4432526_2776    2024_AOAA_Pver_03_1_A_R24196_L002       100.00% 0.00%   0.00%   0.00%
4432526_2777    2024_AOAA_Pver_03_1_A_R24196_L003       4.5     0.0     4.5     0.00%
4432526_2778    2024_AOAA_Pver_03_1_A_R24196_L004       4.5     0.0     4.5     0.00%
4432526_2782    2024_AOAA_Pver_03_1_B_R24196_L008       99.72%  0.28%   0.00%   4.65%
4432526_2779    2024_AOAA_Pver_03_1_B_R24196_L005       99.42%  0.29%   0.29%   6.05%
4432526_2780    2024_AOAA_Pver_03_1_B_R24196_L006       99.72%  0.28%   0.00%   6.13%
4432526_2781    2024_AOAA_Pver_03_1_B_R24196_L007       100.00% 0.00%   0.00%   6.14%
4432526_1500    2024_OFU3_Pver_14_1_A_R24196_L002       50.00%  0.00%   50.00%  50.00%
4432526_1496    2024_OFU3_Pver_13_1_B_R24196_L006       47.57%  33.94%  18.48%  63.32%
4432526_1495    2024_OFU3_Pver_13_1_B_R24196_L005       47.23%  34.17%  18.60%  63.46%
4432526_1498    2024_OFU3_Pver_13_1_B_R24196_L008       47.45%  34.03%  18.53%  63.51%
4432526_1497    2024_OFU3_Pver_13_1_B_R24196_L007       47.20%  34.20%  18.60%  63.55%
4432526_1493    2024_OFU3_Pver_13_1_A_R24196_L003       45.18%  35.34%  19.49%  63.86%
4432526_1494    2024_OFU3_Pver_13_1_A_R24196_L004       45.04%  35.39%  19.57%  63.88%
4432526_1492    2024_OFU3_Pver_13_1_A_R24196_L002       44.93%  35.52%  19.55%  63.99%
4432526_1491    2024_OFU3_Pver_13_1_A_R24196_L001       44.95%  35.50%  19.56%  64.10%
4432526_1160    2024_MALO_Pver_11_1_B_R24196_L006       39.16%  39.10%  21.73%  69.35%
4432526_1159    2024_MALO_Pver_11_1_B_R24196_L005       39.00%  39.20%  21.80%  69.42%
.
.
.
4432526_2753    2024_FASA_Pspp_Extra4_1_A_R24196_L003   20.93%  52.77%  26.30%  89.82%
4432526_2754    2024_FASA_Pspp_Extra4_1_A_R24196_L004   20.78%  52.85%  26.37%  89.84%
4432526_2752    2024_FASA_Pspp_Extra4_1_A_R24196_L002   20.71%  52.93%  26.36%  89.85%
4432526_2751    2024_FASA_Pspp_Extra4_1_A_R24196_L001   20.83%  52.88%  26.28%  89.88%
4432526_3110    2024_FMAL_Pspp_Extra1_1_A_R24196_L004   17.10%  54.17%  28.73%  90.00%
4432526_3109    2024_FMAL_Pspp_Extra1_1_A_R24196_L003   16.67%  54.70%  28.64%  90.23%
4432526_1605    2024_OFU3_Pver_27_1_A_R24196_L003       18.75%  49.43%  31.82%  90.34%
4432526_1502    2024_OFU3_Pver_14_1_A_R24196_L004       16.67%  33.33%  50.00%  91.67%
4432526_1499    2024_OFU3_Pver_14_1_A_R24196_L001       0.00%   100.00% 0.00%   100.00%
4432526_1501    2024_OFU3_Pver_14_1_A_R24196_L003       0.00%   14.29%  85.71%  100.00%
4432526_429     2024_LEON_Pver_11_1_B_R24196_L005       0.00%   100.00% 0.00%   100.00%
4432526_430     2024_LEON_Pver_11_1_B_R24196_L006       0.00%   90.91%  9.09%   100.00%
```

Most libraries had 80-90% mapping rate.
The only libraries with really low mapping were the problematic `AOAA_Pver_03` libraries (sample 823), which had little to no sequence data.
`LEON_Pver_11_1_B` and `OFU3_Pver_14_1_A` also had essentially no data.
`OFU3_Pver_27_1_A_R24196_L003` had only 176 reads.
`2024_FMAL_Pspp_Extra1_1_A` libraries had only about 2000 reads per lane, but B and C libraries looked good.

<br>

Check duplicate rate of each library:
`summarize_dedup.sh`
```bash
#!/bin/bash

# Create output file
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/pver_gwas_batch3/bam/
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

Check highest and lowest duplication rates:
```bash
cut -f1,2,4,10 dupstat_summary.tsv | sort -k4,4 -r | less -S
```

```

FILE    LIBRARY     READ_PAIRS_EXAMINED     PERCENT_DUPLICATION
2024_OFU3_Pver_14_1_A_R24196_L003_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_OFU3_Pver_14       7       0.714286
2024_LEON_Pver_11_1_B_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_LEON_Pver_11       11      0.545455
2024_FTEL_Pver_14_1_B_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_FTEL_Pver_14       10419159        0.435372
2024_OFU6_Pver_21_1_B_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_OFU6_Pver_21       4929755 0.432997
2024_FTEL_Pver_31_1_C_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_FTEL_Pver_31       4901830 0.432441
2024_VATI_Pver_04_1_B_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_VATI_Pver_04       7895108 0.431617
2024_FTEL_Pver_27_1_C_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_FTEL_Pver_27       5560226 0.431447
2024_AOAA_Pver_17_1_B_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_17       7127156 0.43059
2024_FALU_Pver_12_1_B_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_FALU_Pver_12       9529413 0.429456
2024_OFU3_Pver_22_1_B_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_OFU3_Pver_22       6888683 0.428907
.
.
.
2024_LEON_Pver_03_1_B_R24196_L005_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_LEON_Pver_03       1504826 0.135624
2024_LEON_Pver_03_1_B_R24196_L007_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_LEON_Pver_03       1506493 0.134433
2024_LEON_Pver_03_1_B_R24196_L008_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_LEON_Pver_03       1487780 0.132482
2024_AOAA_Pver_03_1_B_R24196_L007_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_03       0       0.046512
2024_OFU3_Pver_14_1_A_R24196_L002_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_OFU3_Pver_14       1       0
2024_LEON_Pver_11_1_B_R24196_L005_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_LEON_Pver_11       1       0
2024_AOAA_Pver_03_1_B_R24196_L008_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_03       1       0
2024_AOAA_Pver_03_1_B_R24196_L006_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_03       1       0
2024_AOAA_Pver_03_1_B_R24196_L005_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_03       2       0
2024_AOAA_Pver_03_1_A_R24196_L004_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_03       0       0
2024_AOAA_Pver_03_1_A_R24196_L003_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_03       0       0
2024_AOAA_Pver_03_1_A_R24196_L002_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_03       0       0
2024_AOAA_Pver_03_1_A_R24196_L001_bt2_combined_pver_cd_hologenome_reheadered_qsorted_dupstat.txt        2024_AOAA_Pver_03       0       0

```
Most samples had **20-40% duplicates**.

<br>



Make list of deduped bam files:
```bash
cd $BASEDIR/pver_gwas/pver_gwas_batch3/bam

# Exclude Ahya samples and AOAA_Pver_03
ls *dedup_coordsorted.bam | grep -v "Ahya" | grep -v "AOAA_Pver_03" > $BASEDIR/pver_gwas/pver_gwas_batch3/sample_lists/first_dedup_bams_list.txt
```
This leaves 3086 BAM files.

<br>

## Create symbolic links (shortcuts) for BAM files in a new directory called `$BASEDIR/pver_gwas/hologenome_mapped_all/`

```bash
BAMDIR=$BASEDIR/pver_gwas/pver_gwas_batch3/bam
BAMLIST=$BASEDIR/pver_gwas/pver_gwas_batch3/sample_lists/first_dedup_bams_list.txt

for FILE in $(cat $BAMLIST); do
  ln -s $BAMDIR/$FILE $BASEDIR/pver_gwas/hologenome_mapped_all/
done
```
There are now a total of 4630 bam files linked to this directory.

<br>

Double check number of unique libraries remaining:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all

ls *.bam | cut -f1-6 -d "_" | sort | uniq | wc -l
```
```
772
```
772 libraries, 14 of which are C libraries.

<br>

Double check number of unique extractions remaining:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all

ls *.bam | cut -f1-5 -d "_" | sort | uniq | wc -l
```
```
380
```

<br>


## Run script `merge_bams_by_sample_array.slurm` to merge bams from the 2 different batches (4 + 8 = 12 lanes per library) by sample
```bash
#!/bin/bash

#SBATCH --job-name=merge_bams_by_sample_2025-06-02
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-380%110
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env samtools

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/hologenome_mapped_all
OUTDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams
SAMPLELIST=($(ls $BAMDIR/*.bam | sed -E 's/.*\/2024_([A-Z]{3}._P[a-z]{3}_[0-9A-Za-z]+_[12])_.*/\1/' | sort | uniq))
REFBASENAME=combined_pver_cd_hologenome


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
crun.samtools samtools merge -f -@ 16 "$MERGEDBAM" $BAMFILES


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


## Query-sort, deduplicate, and then coordinate-sort the merged bams:
`dedup_merged_bams_array.slurm`
```bash
#!/bin/bash

#SBATCH --job-name dedup_merged_bams_array_2025-06-03
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-380%110
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env
module load container_env gatk
GATK='crun.gatk gatk'

## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams
OUTDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams
SAMPLELIST=$BASEDIR/pver_gwas/hologenome_mapped_all/sample_lists/merged_bams_list.txt

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
$GATK --java-options "-Xmx120G" SortSam \
  --INPUT $BAMDIR/$SAMPLEFILE \
  --OUTPUT $BAMDIR/${SAMPLEFILE%.*}'_qsorted.bam' \
  --SORT_ORDER queryname

# Run validation of BAM file
$GATK --java-options "-Xmx120G" ValidateSamFile \
  -I $BAMDIR/${SAMPLEFILE%.*}'_qsorted.bam' \
  -O $BAMDIR/${SAMPLEFILE%.*}'_qsorted.val' \
  -M VERBOSE

## Mark and remove duplicates
$GATK --java-options "-Xmx120G" MarkDuplicates \
  -I $BAMDIR/${SAMPLEFILE%.*}'_qsorted.bam' \
  -O $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.bam' \
  --METRICS_FILE $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dupstat.txt' \
  --REMOVE_DUPLICATES true

## Run validation of BAM file
$GATK --java-options "-Xmx120G" ValidateSamFile \
  -I $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.bam' \
  -O $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.val' \
  -M VERBOSE
  
## Remove qsorted BAM files to make space
rm $BAMDIR/${SAMPLEFILE%.*}'_qsorted.bam'

## Run SortSam to sort by coordinate for downstream processing
$GATK --java-options "-Xmx120G" SortSam \
  --INPUT $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.bam' \
  --OUTPUT $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup_coordsorted.bam' \
  --SORT_ORDER coordinate

## Run validation of coordinate-sorted BAM file
$GATK --java-options "-Xmx120G" ValidateSamFile \
  -I $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup_coordsorted.bam' \
  -O $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup_coordsorted.val' \
  -M VERBOSE
  
## Remove qsorted deduped BAM files to make space
rm $OUTDIR/${SAMPLEFILE%.*}'_qsorted_dedup.bam'
```

<br>

Check duplicate rate of each library:
`summarize_dedup.sh`
```bash
#!/bin/bash

# Create output file
BASEDIR=/archive/barshis/barshislab/jtoy/
BAMDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/
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
**2-5% of reads** were identified as duplicates from each merged bam file.

<br>

### Make list of deduped merged bam files
```bash
cd $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams

ls *qsorted_dedup_coordsorted.bam > ../../sample_lists/dedup_bams_coordsorted_list.txt
```

<br>

## Count alignments remaining post dedup
`count_postdedup_reads_array.slurm`
```bash
#!/bin/bash
#SBATCH --job-name count_postdedup_reads_array_2025-06-04
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-380%110
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env samtools

BASEDIR=/archive/barshis/barshislab/jtoy/
BAMLIST=$BASEDIR/pver_gwas/hologenome_mapped_all/sample_lists/dedup_bams_coordsorted_list.txt

## Change working directory
cd $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams


## Loop over each sample
SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


# Count mapped reads
COUNT=`crun.samtools samtools view -@16 -F 4 $SAMPLEBAM | wc -l`

echo -e "$SAMPLEBAM\t$COUNT" > $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/${SAMPLEBAM%.*}_mappedcount.txt


echo "done-zo woot!"
```


Compile individual counts:
```bash
echo -e "Sample\tMapped_read_count" > postdedup_mapped_read_counts.txt
cat *mappedcount.txt >> postdedup_mapped_read_counts.txt
```

<br>

Check mapped read counts:
```bash
sort -V -k2,2 postdedup_mapped_read_counts.txt | less
```
```
Sample  Mapped_read_count
2024_OFU3_Pver_24_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    49135709
2024_OFU3_Pver_11_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    53602320
2024_OFU3_Pver_05_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    55839027
2024_OFU3_Pver_07_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    58540183
2024_OFU3_Pver_16_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    67105148
2024_OFU3_Pver_12_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    67352324
2024_OFU3_Pver_14_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    71918578
2024_FASA_Pver_34_2_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    77231526
2024_ALOF_Pver_27_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    81733000
2024_ALOF_Pver_26_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    84190748
.
.
.
2024_FALU_Pver_33_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    206807809
2024_AOAA_Pver_14_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    216418108
2024_VATI_Pver_11_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    220963239
2024_VATI_Pver_08_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    228070385
2024_FTEL_Pver_38_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    229796295
2024_OFU3_Pver_27_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    233646488
2024_LEON_Pver_07_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    233696810
2024_MALO_Pver_03_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    233794473
2024_VATI_Pver_15_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    240919552
2024_ALOF_Pver_05_1_combined_pver_cd_hologenome_merged_qsorted_dedup_coordsorted.bam    293207462
```

<br>

## Calculate depth/coverage statistics with CollectWgsMetrics:
`CollectWgsMetrics_array.slurm`
```bash
#!/bin/bash
#SBATCH --job-name CollectWgsMetrics_hologenome_array_2025-06-04
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-380%110
#SBATCH --ntasks=1
#SBATCH --mem=110G
#SBATCH --time 7-00:00:00



## Load modules
module load container_env gatk

BASEDIR=/archive/barshis/barshislab/jtoy/
BAMLIST=$BASEDIR/pver_gwas/hologenome_mapped_all/sample_lists/dedup_bams_coordsorted_list.txt
GATK='crun.gatk gatk'
REFERENCE=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/combined_pver_cd_hologenome.fa


## Loop over each sample

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run CollectWgsMetrics
$GATK --java-options "-Xmx110G" CollectWgsMetrics \
  --INPUT $BASEDIR'/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/'$SAMPLEBAM \
  --OUTPUT $BASEDIR'/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/'${SAMPLEBAM%.*}'_metrics.txt' \
  --REFERENCE_SEQUENCE $REFERENCE

echo 'done-zo woot!'
```

<br>


## Extract only host (Pver) alignments:
Filter bams using `select_host_alignments_array.slurm` script:

```bash
#!/bin/bash

#SBATCH --job-name select_host_alignments_array_2025-06-04
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-380%110
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 10-00:00:00




## Define some variables
BASEDIR=/archive/barshis/barshislab/jtoy/
OUTDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams
SAMPLELIST=$BASEDIR/pver_gwas/hologenome_mapped_all/sample_lists/dedup_bams_coordsorted_list.txt # Path to a bam list
REFBASENAME=PverCD
SCAFLIST=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/Pver_scaffold_names_singleline.txt


# Make new directory
mkdir -p $OUTDIR

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
crun.samtools samtools index -@ 16 $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/$SAMPLEFILE

# Extract mappings that are primary alignments only (no unmapped or secondary/supplementary reads), with mapping score > 20, mapping length (CIGAR) > 20, and only on host scaffolds. Then extract only reads that aligned concordantly (-f 2)
crun.samtools samtools view -b -F 260 -q 20 -m 20 -@ 16 $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/$SAMPLEFILE `cat $SCAFLIST` | crun.samtools samtools view -@ 16 -f 2 -b -o $OUTDIR/$HOSTOUT

# Remove Cgoreaui and Dtrenchii sequence header lines, then reheader bam
crun.samtools samtools view -H $OUTDIR/$HOSTOUT | sed -e '/Cgoreaui/d' -e '/Dtrenchii/d' > $OUTDIR/'header_'$SLURM_ARRAY_TASK_ID'.sam'

crun.samtools samtools reheader $OUTDIR/'header_'$SLURM_ARRAY_TASK_ID'.sam' $OUTDIR/$HOSTOUT > $OUTDIR/${HOSTOUT%.*}'_reheadered.bam'

# Remove original bam file
rm $OUTDIR/$HOSTOUT
```

<br>

Make list of Pver bams:
```bash
cd $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/
ls *.bam > ../../../sample_lists/pver_bams_list.txt
```

<br>

Count host (Pver) reads:
`count_host_reads_array.slurm`
```bash
#!/bin/bash
#SBATCH --job-name count_host_reads_array_2025-06-05
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-380%110
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env samtools

BASEDIR=/archive/barshis/barshislab/jtoy/
BAMLIST=$BASEDIR/pver_gwas/hologenome_mapped_all/sample_lists/pver_bams_list.txt

## Change working directory
cd $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams


## Loop over each sample
SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


# Count mapped reads
COUNT=`crun.samtools samtools view -@16 -F 4 $SAMPLEBAM | wc -l`

echo -e "$SAMPLEBAM\t$COUNT" > $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/${SAMPLEBAM%.*}_pvercount.txt


echo "done-zo woot!"
```

Compile counts
```bash
echo -e "Sample\tPver_mapped_read_count" > pver_mapped_read_counts.txt
cat *pvercount.txt >> pver_mapped_read_counts.txt
```
 <br>
```bash
sort -V -k2,2 pver_mapped_read_counts.txt | less
```

```
Sample  Pver_mapped_read_count
2024_OFU3_Pver_24_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      25308146
2024_OFU3_Pver_05_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      27316694
2024_OFU3_Pver_11_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      28248246
2024_OFU3_Pver_07_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      29778938
2024_OFU3_Pver_16_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      33678558
2024_OFU3_Pver_12_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      34319626
2024_OFU3_Pver_14_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      37190408
2024_FASA_Pver_34_2_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      39548008
2024_ALOF_Pver_27_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      40940074
2024_ALOF_Pver_26_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      42459216
.
.
.
2024_FALU_Pver_33_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      104546770
2024_AOAA_Pver_14_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      113136116
2024_VATI_Pver_11_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      114869950
2024_LEON_Pver_07_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      117219424
2024_FTEL_Pver_38_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      119316852
2024_VATI_Pver_08_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      119600210
2024_OFU3_Pver_27_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      124563044
2024_MALO_Pver_03_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      125612594
2024_VATI_Pver_15_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      130513866
2024_ALOF_Pver_05_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam      154967226
```
 
### Calculate percent host reads
```bash
BASEDIR=/archive/barshis/barshislab/jtoy/

# first sort counts files
tail -n +2 $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/postdedup_mapped_read_counts.txt | sort -k1,1 >  $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/postdedup_mapped_read_counts_sorted.txt
tail -n +2 $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/pver_mapped_read_counts.txt | sort -k1,1 > $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/pver_mapped_read_counts_sorted.txt

# make sure first column of both files are identical
diff <(cut -f1 $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/postdedup_mapped_read_counts_sorted.txt | cut -f1-5 -d "_") \
     <(cut -f1 $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/pver_mapped_read_counts_sorted.txt | cut -f1-5 -d "_")

# paste columns together
paste <(cut -f1 "$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/pver_mapped_read_counts_sorted.txt") \
      <(cut -f2 "$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/pver_mapped_read_counts_sorted.txt") \
      <(cut -f2 "$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/postdedup_mapped_read_counts_sorted.txt") \
      > "$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/pver_and_total_mapped_read_counts.txt"

# calculate percentages
awk '{print $0, $2/$3}' \
    "$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/pver_and_total_mapped_read_counts.txt" \
    > "$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/pver_mapped_read_percents.txt"
```

```bash
sort -k4,4 pver_mapped_read_percents.txt | less
```
| Sample                                                                         | Pver_mapped_read_count | Total_mapped_read_count | Percent_pver_reads |
|--------------------------------------------------------------------------------|------------------------|-------------------------|--------------------|
| 2024_LEON_Pver_11_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam     |           51,437,792   |             106,691,465 | 0.482117           |
| 2024_OFU3_Pver_05_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam     |           27,316,694   |              55,839,027 | 0.489204           |
| 2024_ALOF_Pspp_Extra1_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam |           62,531,432   |             127,636,414 | 0.489918           |
| 2024_AOAA_Pspp_Extra7_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam |           48,211,304   |              98,394,731 | 0.489979           |
| 2024_FTEL_Pver_13_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam     |           53,235,276   |             108,582,167 | 0.490276           |
| 2024_AOAA_Pver_34_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam     |           68,500,054   |             139,624,415 | 0.490602           |
| 2024_MALO_Pver_26_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam     |           55,815,728   |             113,693,234 | 0.490933           |
| 2024_AOAA_Pver_18_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam     |           44,944,396   |              91,530,559 | 0.491032           |
| 2024_FTEL_Pspp_Extra4_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam |           84,261,370   |             171,551,884 | 0.491171           |
| 2024_AOAA_Pver_17_1_PverCD_dedup_primary_minq20_mlen20_pver_reheadered.bam     |           51,694,440   |             105,192,387 | 0.491428           |

Percent of mapped reads that concordantly map to Pver contigs ranges from **48-60%**.

<br>


## Calculate depth/coverage statistics for host with `CollectWgsMetricsWithNonZeroCoverage`
`CollectWgsMetrics_pver_only_nonzero_array.slurm`
```bash
#!/bin/bash
#SBATCH --job-name CollectWgsMetrics_pver_only_nonzero_array_2025-06-05
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-380%110
#SBATCH --ntasks=1
#SBATCH --mem=110G
#SBATCH --time 7-00:00:00


## Load modules
module load container_env gatk

BASEDIR=/archive/barshis/barshislab/jtoy/
BAMLIST=$BASEDIR/pver_gwas/hologenome_mapped_all/sample_lists/pver_bams_list.txt
GATK='crun.gatk gatk'
REFERENCE=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta


## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Run CollectWgsMetrics
$GATK --java-options "-Xmx110G" CollectWgsMetricsWithNonZeroCoverage \
  --INPUT $BASEDIR'/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/'$SAMPLEBAM \
  --OUTPUT $BASEDIR'/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/'${SAMPLEBAM%.*}'_metrics_nonzero.txt' \
  --CHART $BASEDIR'/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/'${SAMPLEBAM%.*}'_metrics_nonzero_chart.pdf' \
  --REFERENCE_SEQUENCE $REFERENCE

echo 'done-zo woot!'
```

<br>

### Plot histograms of coverage depth and summarise coverage statistics using R script `summarize_coverage_stats.R`
First extract histogram data from output files and save to new files:
```bash
for FILE in `ls *nonzero.txt`; do
    tail -n +12 $FILE > ${FILE%.*}_histdata.tsv
done
```

<br>

I realized some files have "Extra" in their name which throws off the consistency of file name character lengths. So first rename these files, replacing "Extra" with "X" so that this part of the sample identifier becomes 2 characters in length.
```bash
# First do a dry run to make sure the command is working as expected:
for FILE in *Extra*_PverCD_dedup_primary_minq20_mlen20_pver_reheadered*; do
  newname=$(echo "$FILE" | sed 's/Extra/X/')
  echo $newname
done

# Now actually change the names
for FILE in *Extra*_PverCD_dedup_primary_minq20_mlen20_pver_reheadered*; do
  newname=$(echo "$FILE" | sed 's/Extra/X/')
  mv "$FILE" "$newname"
done
```

Then run R script (I ran it in Rstudio Server to track it in real time):
`summarize_coverage_stats.R`
```r
setwd("/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/")

file_list <- list.files(pattern = "histdata\\.tsv$")

for (FILE in file_list) {
  
  # extract sample name
  sample <- substr(FILE, 1, 19)
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

Check per-sample coverage/depth:
```bash
sort -k3,3 coverage_summary.tsv | less
```

| Sample              | Mean_Depth                                               | Median_Depth                          | Proportion_Genome_Covered |
|---------------------|----------------------------------------------------------|---------------------------------------|---------------------------|
| 2024_OFU3_Pver_24_1 |                                                  13.806  |                                   14  | 0.670632                  |
| 2024_OFU3_Pver_11_1 |                                                  14.312  |                                   15  | 0.672853                  |
| 2024_OFU3_Pver_05_1 |                                                  15.489  |                                   16  | 0.675293                  |
| 2024_OFU3_Pver_07_1 |                                                  16.627  |                                   17  | 0.678445                  |
| 2024_OFU3_Pver_12_1 |                                                  18.469  |                                   19  | 0.682938                  |
| 2024_OFU3_Pver_16_1 |                                                  18.712  |                                   19  | 0.68163                   |
| 2024_OFU3_Pver_14_1 |                                                  20.102  |                                   21  | 0.683576                  |
| 2024_ALOF_Pver_27_1 |                                                  22.375  |                                   23  | 0.68571                   |
| 2024_FASA_Pver_34_2 |                                                  21.893  |                                   23  | 0.686843                  |
| 2024_AOAA_Pver_02_1 |                                                  23.234  |                                   24  | 0.687666                  |
|...                  |...                                                       |...           |...                        |
| 2024_VATI_Pver_11_1 |                                                  53.026  | 59           | 0.710328                  |
| 2024_FALU_Pver_33_1 |                                                  54.787  | 61           | 0.708923                  |
| 2024_LEON_Pver_07_1 |                                                  57.505  | 64           | 0.711158                  |
| 2024_AOAA_Pver_14_1 |                                                  57.742  | 65           | 0.711169                  |
| 2024_FTEL_Pver_38_1 |                                                  59.097  | 67           | 0.744652                  |
| 2024_VATI_Pver_08_1 |                                                  60.125  | 67           | 0.715447                  |
| 2024_MALO_Pver_03_1 |                                                  60.232  | 68           | 0.704                     |
| 2024_OFU3_Pver_27_1 |                                                  62.156  | 70           | 0.713505                  |
| 2024_VATI_Pver_15_1 |                                                  65.637  | 74           | 0.710618                  |
| 2024_ALOF_Pver_05_1 |                                                  76.069  | 87           | 0.716447                  |

Most samples have at least **20x median** depth and **~70% coverage**.

<br>

### Create new list of bam files that incorporates the name changes
```bash
cd $BASEDIR/jtoy/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/

ls *.bam > ../../../sample_lists/pver_bams_list_renamed.txt
```

<br>

## Call variants per sample with HaplotypeCaller (in GVCF mode)

To run HaplotypeCaller, the reference fasta must first be indexed (samtools) and a sequence dictionary created (GATK/picard). This was already done for the analysis of the pilot data using the commands below:
```bash
# samtools indexing
module load samtools
cd /cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1
crun.samtools samtools faidx GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta

# create sequence dictionary with GATK/picard
module load gatk
GATK='crun.gatk gatk'
$GATK --java-options "-Xmx100G" CreateSequenceDictionary --REFERENCE GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta
```

Run HaplotypeCaller with the following array script:
`HaplotypeCaller_pver_array`
```bash
#!/bin/bash
#SBATCH --job-name HaplotypeCaller_pver_array_2025-06-06
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-380%110
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 14-00:00:00
#SBATCH --cpus-per-task=16


## Load modules
module load container_env gatk

BASEDIR=/archive/barshis/barshislab/jtoy/
BAMLIST=$BASEDIR/pver_gwas/hologenome_mapped_all/sample_lists/pver_bams_list_renamed.txt
GATK='crun.gatk gatk'
REFERENCE=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF>OUTDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/gvcfs/

## Get sample BAM filename
SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BAMLIST)
BAMFILE=$BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/$SAMPLEBAM
BAIFILE=${BAMFILE%.*}.bai

echo "Slurm array task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample BAM: $SAMPLEBAM"

## Check if BAM index exists, and create it if missing
if [[ ! -f "$BAIFILE" ]]; then
  echo "Index file $BAIFILE not found. Creating index..."
  cd $BASEDIR/pver_gwas/hologenome_mapped_all/merged_bams/dedup_bams/pver_bams/
  $GATK --java-options "-Xmx25G" BuildBamIndex --INPUT "$BAMFILE"
else
  echo "Index file $BAIFILE already exists. Skipping indexing."
fi

## Make directory for GVCFs if it doesn't already exist
mkdir -p $OUTDIR

## Run HaplotypeCaller
$GATK --java-options "-Xmx25G" HaplotypeCaller \
  -I $BAMFILE \
  -O $OUTDIR/${SAMPLEBAM%.*}'.g.vcf.gz' \
  -R $REFERENCE \
  -ERC GVCF \
  --native-pair-hmm-threads 16

echo "done-zo woot!"
```
The current max per-user CPU usage for Wahab is 512, so using 16 threads per job allows 32 jobs to run simultaneously. Each job seems to take anywhere from 11-20 hrs to complete with most around 13. Total run time for the array job was 7 days.

<br>

Checked CPU use efficiency using:
```bash
sacct -j 4457047 --format=JobID,JobName%25,AllocCPUs,Elapsed,TotalCPU,CPUTimeRAW,MaxRSS,State
```
Calculated efficency as: TotalCPU / (AllocCPUS * Elapsed). Efficiency for jobs was around 50-60%.

<br>

## Consolidate per sample GVCFs with GenomicsDBImport
First need to create a tab-delimited map file with sample_name in the first column and /path/to/vcf in the second column:
```bash
BASEDIR=/archive/barshis/barshislab/jtoy/
GVCFDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/gvcfs/

cd $GVCFDIR

OUTPUT="pver_gwas_batch23_gvcf.sample_map"
> $OUTPUT

# Loop through matching files
for FILE in `ls $GVCFDIR*_reheadered.g.vcf.gz`; do
  # Extract first 5 underscore-separated fields
  PREFIX=$(echo "$FILE" | cut -d '_' -f1-5)
  
  # Append to the output file
  echo -e "$PREFIX\t$GVCFDIR$FILE" >> $OUTPUT
done
```

Then combine this map file with the map file for the pilot samples:
```bash
cat pver_gwas_batch23_gvcf.sample_map /cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/gvcfs/pver_gwas_pilot_gvcf.sample_map | sort -k1,1 > pver_gwas_all_gvcf.sample_map
```
This map file should now have 396 entries.


### Run GenomicsDBImport SLURM script
`GenomicsDBImport.slurm`:
```bash
#!/bin/bash

#SBATCH --job-name GenomicsDBImport_pver_2025-06-17
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --mem=340G
#SBATCH --time 7-00:00:00
#SBATCH --cpus-per-task=38

## Load modules
module load container_env gatk

BASEDIR=/archive/barshis/barshislab/jtoy/
GATK='crun.gatk gatk'

## Run GenomicsDBImport
$GATK --java-options "-Xmx300g -Xms10g" \
   GenomicsDBImport \
   --genomicsdb-workspace-path $BASEDIR/pver_gwas/hologenome_mapped_all/genomicsdb/ \  # Note: do not create this directory ahead of time. GATK will create it automatically
   --sample-name-map $BASEDIR/pver_gwas/hologenome_mapped_all/gvcfs/pver_gwas_all_gvcf.sample_map \
   -L /cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/genome_regions.list \
```
This was too slow, so I canceled it and rescripted it as an array job, splitting by genomic region (contig).

`GenomicsDBImport_array.slurm`
```bash
#!/bin/bash

#SBATCH --job-name GenomicsDBImport_array_pver_2025-06-18
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 7-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --array=1-52%51

## Load modules
module load container_env gatk
GATK='crun.gatk gatk'


## Assign variables
BASEDIR=/archive/barshis/barshislab/jtoy/
REGION_LIST=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_0366699>

## Create genomicsdb directory if it doesn't already exist
mkdir -p $BASEDIR/pver_gwas/hologenome_mapped_all/genomicsdb


## Initiate array
REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $REGION_LIST)


## Run GenomicsDBImport
$GATK --java-options "-Xmx90g -Xms10g" \
   GenomicsDBImport \
   --genomicsdb-workspace-path $BASEDIR/pver_gwas/hologenome_mapped_all/genomicsdb/region_${SLURM_ARRAY_TASK_ID} \
   --sample-name-map $BASEDIR/pver_gwas/hologenome_mapped_all/gvcfs/pver_gwas_all_gvcf.sample_map \
   -L $REGION \
   --reader-threads 10
```
The first chromosome took the longest to import - a little over 24 hrs.

<br>

## Joint genotyping with GenotypeGVCFs
This step is also parallelized by genomic region to speed up processing time, since GenotypeGVCFs doesn't have internal parallelization options.

`GenotypeGVCFs_array.slurm`
```bash
#!/bin/bash

#SBATCH --job-name=GenotypeGVCFs_pver_2025-06-19
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-52%51
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=10

# Load GATK module
module load container_env gatk

# Define paths
BASEDIR=/archive/barshis/barshislab/jtoy
REFERENCE=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta
OUTDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/vcf
GENDBBASE=$BASEDIR/pver_gwas/hologenome_mapped_all/genomicsdb
SCAFLIST=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/genome_regions.list
GATK='crun.gatk gatk'

# Create output directory if it doesn't exist
mkdir -p $OUTDIR

# Get the current region and scaffold
REGION_ID=${SLURM_ARRAY_TASK_ID}
GENDB=${GENDBBASE}/region_${REGION_ID}
SCAF=$(sed -n "${REGION_ID}p" $SCAFLIST)

echo "Processing region: region_${REGION_ID}, corresponding to scaffold: $SCAF"
echo "Using GenomicsDB workspace: $GENDB"

# Run GenotypeGVCFs
$GATK --java-options "-Xmx95g" GenotypeGVCFs \
  -R $REFERENCE \
  -V gendb://$GENDB \
  -O $OUTDIR/${SCAF}_genotypes.vcf.gz
```

Checking the .err files, there are a few warning messages that appear:
```
Chromosome NC_089321.1_Pverrucosa position 2941072 (TileDB column 251074565) has too many alleles in the combined VCF record : 62 : current limit : 50. Fields, such as  PL, with length equal to the number of genotypes will NOT be added for this location.
02:00:43.651 WARN  MinimalGenotypingEngine - Attempting to genotype more than 50 alleles. Site will be skipped at location NC_089321.1_Pverrucosa:2941072

Sample/Callset 2024_ALOF_Pspp_X1_1( TileDB row idx 0) at Chromosome NC_089321.1_Pverrucosa position 3250313 (TileDB column 251383806) has too many genotypes in the combined VCF record : 1176 : current limit : 1024 (num_alleles, ploidy) = (48, 2). Fields, such as  PL, with length equal to the number of genotypes will NOT be added for this sample for this location.

02:15:04.517 WARN  MinimalGenotypingEngine - No genotype contained sufficient data to recalculate site and allele qualities. Site will be skipped at location NC_089321.1_Pverrucosa:3250313
```
These are all related and expected warnings for WGS data in a nonmodel organism. They are caused by hyper-variable sites with more than 50 different alleles, likely caused by issues with the reference sequence, misalignments, or repeat regions. These likely low-quality and multiallelic sites will get filtered out of the VCF downsteam anyway. These warnings occurred 318, 101376, and 256 times respectively across all .err files for the array job. 


Once genotyping is complete, there will we a separate vcf for each genomic region (scaffold/contig). Combine them into a single multi-sample VCF with BCFtools:
```bash
# make list of vcf files
ls *.vcf.gz | sort > vcf_list.txt


# combine VCFs
crun.bcftools bcftools concat -f vcf_list.txt --threads 34 -Ov -o pver_all_combined_genotypes.vcf
```

Now count the total number of called genotypes (variants):
```bash
grep -v "#" pver_all_combined_genotypes.vcf | wc -l
```
**32,498,469** total SNPs

<br>

## 12. Variant filtering

First we need to know a bit about the depth of coverage at each variant
```bash
# Get the total site depth per position, i.e., the sum of read depths from all samples at each variant site, using the INFO/DP field.
crun.bcftools bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' pver_all_combined_genotypes.vcf > pver_all_combined_genotypes.total_site_depth                # run on full vcf
```
```
NC_089312.1_Pverrucosa  3708    66
NC_089312.1_Pverrucosa  3822    6
NC_089312.1_Pverrucosa  5798    23
NC_089312.1_Pverrucosa  5823    8
NC_089312.1_Pverrucosa  6389    2
NC_089312.1_Pverrucosa  7705    150
NC_089312.1_Pverrucosa  7712    155
NC_089312.1_Pverrucosa  7719    129
NC_089312.1_Pverrucosa  10033   6
NC_089312.1_Pverrucosa  10799   9
...
```

```bash
# Get the per sample depth at each position
crun.bcftools bcftools query -f '%CHROM\t%POS[\t%DP]\n' pver_all_combined_genotypes.vcf > pver_all_combined_genotypes.per_sample_depth                # run on full vcf
```
```
NC_089312.1_Pverrucosa  3708    0       0       0       0       0       0       0       0       0       2       0       0       0   
NC_089312.1_Pverrucosa  3822    0       0       0       0       0       0       0       0       0       0       0       0       0   
NC_089312.1_Pverrucosa  5798    0       0       0       0       0       0       0       1       0       0       0       0       0   
NC_089312.1_Pverrucosa  5823    0       0       0       0       0       0       0       1       0       0       0       0       0   
NC_089312.1_Pverrucosa  6389    0       0       0       0       0       0       0       0       0       0       0       0       0   
NC_089312.1_Pverrucosa  7705    2       0       0       0       0       0       0       0       0       0       2       0       0   
NC_089312.1_Pverrucosa  7712    2       0       0       0       0       0       0       0       0       0       2       0       0   
NC_089312.1_Pverrucosa  7719    2       0       0       0       0       0       0       0       0       0       3       0       0   
NC_089312.1_Pverrucosa  10033   0       0       0       0       0       0       0       0       0       0       0       0       0   
NC_089312.1_Pverrucosa  10799   0       1       0       0       0       0       0       0       0       0       0       0       0
...
```

```
# Get a set of summary stats for the vcf (distributions of quality scores, depths, etc)
crun.bcftools bcftools stats pver_all_combined_genotypes.vcf > pver_all_combined_genotypes.stats                # run on full vcf
```

<br>

Plot depth per site in R:
```r
# SNP depth of coverage analysis - Pver All Samples
# 2025-07-01
# Jason A. Toy

library(tidyverse)

rm(list = ls())
setwd("/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf")


# load total depth file
td <- read_tsv("pver_all_combined_genotypes.total_site_depth", col_names = FALSE)

# calculate summary stats
mean_td <- mean(td$X3)
median_td <- median(td$X3)
sd_td <- sd(td$X3)


# plot distribution
p <- ggplot(td) +
  geom_freqpoly(aes(x = X3), binwidth = 1) + 
  xlab("SNP depth") +
  xlim(c(0,60000)) +
  ylab("Count") +
  geom_vline(xintercept = mean_td, color = "blue") +
  geom_vline(xintercept = median_td, color = "darkgreen") +
  geom_vline(xintercept = median_td + sd_td, color = "black", linetype = "dashed") +
  annotate("text", x=40000, y=12000, label=paste0(round(length(td$X3)/10^6,2), " million SNPs before filtering"), color = 'red') +
  annotate("text", x=40000, y=11500, label=paste0("mean total depth = ", round(mean_td, 1)), color = "blue") +
  annotate("text", x=40000, y=11000, label=paste0("median total depth = ", median_td), color = "darkgreen") +
  annotate("text", x=40000, y=10500, label=paste0("--- ", "median total depth + 1 SD = ", round(median_td + sd_td, 1)), color = "black") +
  theme_bw()

ggsave(p, "snp_depth_plot_pver_all.png", width = "12", height = "8", units = "in")
```
![](https://github.com/user-attachments/assets/45b1c578-305e-4af1-8440-399d1d9023e6)

<br>

### Filter variants with with PLINK2
First filter on quality scores (mapping & variant) and depth with BCFtools:
```bash
crun.bcftools bcftools filter --threads 36 -e 'QUAL < 30 || INFO/MQ < 40 || INFO/DP < 792 || INFO/DP > 25286' pver_all_combined_genotypes.vcf -Oz -o pver_all_QDPfiltered_genotypes.vcf.gz

# Base quality score filter set to 30. Mapping quality score filter set to 40
# There are 396 total samples now. Per sampl depth filter was set to 2 * 396 = 792. Total depth filter was set to median total depth + 1 SD = 25,286.
```

Count remaining SNPs:
```bash
zcat pver_all_QDPfiltered_genotypes.vcf.gz | grep -v "#" | wc -l
```
This leaves **5,255,644** SNPs

<br>

Next filter based on missingess and MAF and remove indels and multiallelic SNPs:
```bash
module load plink/2024.03.02    # Note: this is PLINK2

crun.plink plink2 \
  --vcf pver_all_QDPfiltered_genotypes.vcf.gz \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno 0.2 \
  --mind 0.2 \
  --maf 0.05 \
  --make-pgen \
  --out pver_all_MISSMAFfiltered_genotypes


  # geno 0.2    removes SNPs where more than 20% of individuals have missing genotypes
  # mind 0.2    removes individuals who are missing more than 20% of genotype data across all SNPs
  # maf 0.05    removes SNPs where the minor allele frequency is less than 0.05
  # threads N    multithreading option available
```
```
...
0 samples removed due to missing genotype data (--mind).
396 samples (0 females, 0 males, 396 ambiguous; 396 founders) remaining aftermain filters.
Calculating allele frequencies... done.
--geno: 8 variants removed due to missing genotype data.
3451112 variants removed due to allele frequency threshold(s)
```

Convert to VCF for viewing:
```bash
crun.plink plink2 \
  --pgen pver_all_MISSMAFfiltered_genotypes.pgen \
  --psam pver_all_MISSMAFfiltered_genotypes.psam \
  --pvar pver_all_MISSMAFfiltered_genotypes.pvar \
  --recode vcf \
  --out pver_all_MISSMAFfiltered_genotypes
```
This leaves **1,134,549** SNPs remaining after filtering.


## LD-pruning with PLINK2
PLINK2 requires unique values in the ID field for each SNP. Our dataset does not contain any ID values, but instead has "." placeholders. So the first thing we need to do is replace these with a chromosome/position-based ID:
```bash
crun.plink plink2 \
  --pgen pver_all_MISSMAFfiltered_genotypes.pgen \
  --psam pver_all_MISSMAFfiltered_genotypes.psam \
  --pvar pver_all_MISSMAFfiltered_genotypes.pvar \
  --set-all-var-ids @:#:\$r:\$a \
  --make-pgen \
  --out pver_all_MISSMAFfiltered_uniqIDs
```
The `@:#:$r:$a` syntax creates IDs formatted as `CHROM:POS:REF:ALT`.

<br>

### Investigate rate of LD decay
Compute pairwise-LD (r2) values:
```
crun.plink plink2 \
  --pfile pver_all_MISSMAFfiltered_uniqIDs \
  --thin-count 100000 \
  --r2-unphased \
  --ld-window-kb 1000 \
  --ld-window-r2 0 \
  --threads 30 \
  --out ld_decay
```

Plot LD decay in R:
`plot_ld_decay.R`
```r
# Plot LD decay - Pver all samples
# 2025-07-02
# Jason A. Toy

library(tidyverse)

rm(list = ls())
setwd("/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf")


# Import LD values
library(data.table) # for faster import of large datasets
ld <- fread("ld_decay.vcor") %>% mutate(DIST_BP = abs(POS_B - POS_A))


# Bin distances and summarize r2 by bin
bin_size <- 1000  # define distance bins in bp

ld_summary <- ld %>%
  mutate(DIST_BIN = floor(DIST_BP / bin_size) * bin_size) %>%
  group_by(DIST_BIN) %>%
  summarise(mean_r2 = mean(UNPHASED_R2, na.rm = TRUE), .groups = "drop")


# Plot LD decay
ldplot <- ggplot(ld_summary, aes(x = DIST_BIN/1000, y = mean_r2)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 0.5, alpha = 0.6) +
  labs(x = "Distance between SNPs (kb)", y = expression(paste("Mean ", r^2)),
       title = "LD Decay Curve - Pver") +
  theme_minimal(base_size = 14)

ggsave(ldplot, file = "ld_decay_plot.pdf")
```
![](https://github.com/user-attachments/assets/75bbc400-7b23-4e7b-b3fc-3fff19f10c57)
The plot shows that **LD declines well below 0.2 by 10 kb**, so I will use that as my window size for pruning.

<br>

### Run LD pruning
Output is an LD-pruned list of SNPs. I ran this twice with different r2 values for comparison:

Stringent pruning for PCA, ADMIXTURE, etc (pruning SNP pairs with r2 > 0.2):
```bash
crun.plink plink2 \
  --pgen pver_all_MISSMAFfiltered_uniqIDs.pgen \
  --psam pver_all_MISSMAFfiltered_uniqIDs.psam \
  --pvar pver_all_MISSMAFfiltered_uniqIDs.pvar \
  --indep-pairwise 10kb 1 0.2 \
  --out pver_all_MISSMAFfiltered_ld_0.2
```
**80,208** SNPs retained.

Lenient pruning (pruning SNP pairs with r2 > 0.5; balancing redundancy minimization with signal retention):
```bash
crun.plink plink2 \
  --pgen pver_all_MISSMAFfiltered_uniqIDs.pgen \
  --psam pver_all_MISSMAFfiltered_uniqIDs.psam \
  --pvar pver_all_MISSMAFfiltered_uniqIDs.pvar \
  --indep-pairwise 10kb 1 0.5 \
  --out pver_all_MISSMAFfiltered_ld_0.5
```
**189,630** SNPs retained

`indep-pairwise` parameters used:
- `10kb` → window size in kb (can also be specified in number of SNPs by omitting "kb" suffix)
- `1` → step size (how many SNPs to shift the window each time; must be 1 bp when specifying window size in kb)
- `0.2` → r² threshold (SNPs with pairwise r² > 0.2 are considered linked)

I used r2 > 0.2 to conservatively remove linked SNPs for population structure analysis, but the more lenient SNP set (r2 > 0.5) may be useful for other downstream analyses.

### Filter dataset using pruned SNP list (stringent)
```bash
crun.plink plink2 \
  --pgen pver_all_MISSMAFfiltered_uniqIDs.pgen \
  --pvar pver_all_MISSMAFfiltered_uniqIDs.pvar \
  --psam pver_all_MISSMAFfiltered_uniqIDs.psam \
  --extract pver_all_MISSMAFfiltered_ld_0.2.prune.in \
  --make-pgen \
  --out pver_all_ld_pruned_0.2_genotypes
```

Convert to VCF for viewing/downstream use:
```bash
crun.plink plink2 \
  --pgen pver_all_ld_pruned_0.2_genotypes.pgen \
  --psam pver_all_ld_pruned_0.2_genotypes.psam \
  --pvar pver_all_ld_pruned_0.2_genotypes.pvar \
  --recode vcf \
  --out pver_all_ld_pruned_0.2_genotypes
```

Double-check number of SNPs in VCF:
```bash
grep -v "^#" pver_all_ld_pruned_0.2_genotypes.vcf | wc -l
```
```
80208
```


## Principle components analysis

Use PLINK2 to run a PCA:
```bash
crun.plink plink2 \
  --pgen pver_all_ld_pruned_0.2_genotypes.pgen \
  --psam pver_all_ld_pruned_0.2_genotypes.psam \
  --pvar pver_all_ld_pruned_0.2_genotypes.pvar \
  --pca 10 \
  --out pver_all_ld_pruned_0.2_pca
```
`--pca 10` calculates the first 10 principal components using exact PCA. For fast approximate PCA (very efficient, but less accurate), use pca approx 10.

Plot PCA in R:
`plot_PCA_pver_all_samples.R`
```r
# Plot PCAs from VCFs - Pver all samples
# 2025-07-02
# Jason A. Toy

rm(list = ls())

setwd("/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf")

library(dplyr)
library(ggplot2)
library(ggrepel)

# Load eigenvectors and eigenvalues
eigenvec <- read.delim("pver_all_ld_pruned_0.2_pca.eigenvec", header = TRUE, sep = "\t")
eigenval <- read.delim("pver_all_ld_pruned_0.2_pca.eigenval", header = FALSE)

# Calculate percentage of variance explained
percent_var <- (eigenval$V1 / sum(eigenval$V1)) * 100
# 15.831364 12.598753 11.226283 10.414841  9.833342  9.137878  8.641080  8.411272  7.321766  6.583421

# Scree Plot
scree <- percent_var %>%
  as_tibble %>%
  rownames_to_column() %>%
  dplyr::rename(PC = rowname) %>% 
  mutate(PC = as.integer(PC)) %>% 
  mutate(cumulative = cumsum(value))

ggplot(scree) +
  geom_point(aes(x = PC, y = value, color = "per PC")) + 
  geom_line(aes(x = PC, y = value, color = "per PC")) + 
  geom_point(aes(x = PC, y = cumulative/10, color = "cumulative/10"), size = 1) +
  geom_line(aes(x = PC, y = cumulative/10, color = "cumulative/10"), linetype = "dashed") +
  scale_color_manual(name = "", values = c("per PC" = "darkred", "cumulative/10" = "lightblue")) +
  ylab("Percent variance explained") +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL) +
  scale_y_continuous(breaks = 1:16, minor_breaks = NULL) +
  theme_minimal()


# Modify eigenvec for plotting
eigenvec_plot <- eigenvec %>%
  separate(X.IID, into = c("Year", "Location", "Species", "Genotype", "Lib_ID"), sep = "_", extra = "merge") %>% 
  mutate(Geno_ID = paste0(Location, "_", Species, "_", Genotype)) %>% 
  mutate(Species = as.factor(Species),
         Location = as.factor(Location))


# Plot PCAs
library(ggsci)
mypal <- pal_d3("category20")(11)

# PC 1/2
PC12 <- ggplot(eigenvec_plot, aes(x = PC1, y = PC2, shape = Location, color = Species)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.3) +
  scale_shape_manual(values = c(0:10)) +
  xlab(paste0("PC1: ", round(percent_var[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2], 2), "% variance")) +
  #geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  theme_bw()
PC12

# color by species, all dots
PC12_noshape <- ggplot(eigenvec_plot, aes(x = PC1, y = PC2, color = Species)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 1.5, alpha = 0.3) +
  xlab(paste0("PC1: ", round(percent_var[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2], 2), "% variance")) +
  theme_bw()
PC12_noshape

# color by location, all dots, shape by species
PC12_color <- ggplot(eigenvec_plot, aes(x = PC1, y = PC2, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  scale_shape_manual(values = c(15, 16)) +
  geom_point(size = 2, alpha = 0.3) +
  xlab(paste0("PC1: ", round(percent_var[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2], 2), "% variance")) +
  #geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  theme_bw()
PC12_color

# faceted by species
PC12_species <- ggplot(eigenvec_plot, aes(x = PC1, y = PC2, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 1.5, alpha = 0.3) +
  xlab(paste0("PC1: ", round(percent_var[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2], 2), "% variance")) +
  #geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Species) +
  theme_bw()

# faceted by location
PC12_location <- ggplot(eigenvec_plot, aes(x = PC1, y = PC2, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 1.5, alpha = 0.3) +
  xlab(paste0("PC1: ", round(percent_var[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2], 2), "% variance")) +
  #geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Location) +
  theme_bw()
PC12_location




# PC 3/4
PC34 <- ggplot(eigenvec_plot, aes(x = PC3, y = PC4, shape = Location, color = Species)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 1.5, alpha = 0.3) +
  scale_shape_manual(values = c(0:10)) +
  xlab(paste0("PC3: ", round(percent_var[3], 2), "% variance")) +
  ylab(paste0("PC4: ", round(percent_var[4], 2), "% variance")) +
  #geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  theme_minimal()
PC34

# color by species, all dots
PC34_noshape <- ggplot(eigenvec_plot, aes(x = PC3, y = PC4, color = Species)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 1.5, alpha = 0.3) +
  xlab(paste0("PC3: ", round(percent_var[3], 2), "% variance")) +
  ylab(paste0("PC4: ", round(percent_var[4], 2), "% variance")) +
  #geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  theme_minimal()
PC34_noshape

# color by location, all dots, shape by species
PC34_color <- ggplot(eigenvec_plot, aes(x = PC3, y = PC4, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  scale_shape_manual(values = c(15, 16)) +
  geom_point(size = 2, alpha = 0.3) +
  xlab(paste0("PC3: ", round(percent_var[3], 2), "% variance")) +
  ylab(paste0("PC4: ", round(percent_var[4], 2), "% variance")) +
  #geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  theme_bw()
PC34_color

# faceted by species
PC34_species <- ggplot(eigenvec_plot, aes(x = PC3, y = PC4, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC3: ", round(percent_var[3], 2), "% variance")) +
  ylab(paste0("PC4: ", round(percent_var[4], 2), "% variance")) +
  #geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Species) +
  theme_bw()
PC34_species

# faceted by location
PC34_location <- ggplot(eigenvec_plot, aes(x = PC3, y = PC4, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC3: ", round(percent_var[3], 2), "% variance")) +
  ylab(paste0("PC4: ", round(percent_var[4], 2), "% variance")) +
  geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Location) +
  theme_bw()
PC34_location




# PC56
PC56 <- ggplot(eigenvec_plot, aes(x = PC5, y = PC6, shape = Location, color = Species)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 1.5, alpha = 0.3) +
  scale_shape_manual(values = c(0:10)) +
  xlab(paste0("PC5: ", round(percent_var[5], 2), "% variance")) +
  ylab(paste0("PC6: ", round(percent_var[6], 2), "% variance")) +
  #geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  theme_minimal()
PC56

# color by species, all dots
PC56_noshape <- ggplot(eigenvec_plot, aes(x = PC5, y = PC6, color = Species)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 1.5, alpha = 0.3) +
  xlab(paste0("PC5: ", round(percent_var[5], 2), "% variance")) +
  ylab(paste0("PC6: ", round(percent_var[6], 2), "% variance")) +
  #geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  theme_minimal()
PC56_noshape

# color by location, all dots, shape by species
PC56_color <- ggplot(eigenvec_plot, aes(x = PC5, y = PC6, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  scale_shape_manual(values = c(15, 16)) +
  geom_point(size = 2, alpha = 0.3) +
  xlab(paste0("PC5: ", round(percent_var[5], 2), "% variance")) +
  ylab(paste0("PC6: ", round(percent_var[6], 2), "% variance")) +
  #geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  theme_bw()
PC56_color

# faceted by species
PC56_species <- ggplot(eigenvec_plot, aes(x = PC5, y = PC6, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC5: ", round(percent_var[5], 2), "% variance")) +
  ylab(paste0("PC6: ", round(percent_var[6], 2), "% variance")) +
  #geom_text_repel(aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Species) +
  theme_bw()
PC56_species

# faceted by location
PC56_location <- ggplot(eigenvec_plot, aes(x = PC5, y = PC6, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC5: ", round(percent_var[5], 2), "% variance")) +
  ylab(paste0("PC6: ", round(percent_var[6], 2), "% variance")) +
  #geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Location) +
  theme_bw()
PC56_location





# combo plots
library(cowplot)
plot_grid(PC12_location, PC34_location, ncol = 2)
plot_grid(PC12_species, PC34_species, ncol = 2)
plot_grid(PC12, PC34, ncol = 2)
plot_grid(PC12_noshape, PC34_noshape)
plot_grid(PC12_color, PC34_color, PC56_color, ncol = 3)





# 3D plots
library(plotly)

# PCs 1-3
plot_ly(
  data = eigenvec_plot,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Location,
  colors = mypal,
  symbol = ~Species,
  symbols = c("square", "circle"),
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


# PCs 4-6
plot_ly(
  data = eigenvec_plot,
  x = ~PC4,
  y = ~PC5,
  z = ~PC6,
  color = ~Location,
  colors = mypal,
  symbol = ~Species,
  symbols = c("square", "circle"),
  text = ~Geno_ID,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 3)
) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("PC4: ", round(percent_var[4], 2), "% variance")),
      yaxis = list(title = paste0("PC5: ", round(percent_var[5], 2), "% variance")),
      zaxis = list(title = paste0("PC6: ", round(percent_var[6], 2), "% variance"))
    )
  )





# Plot various PC combinations
# PC 3/5
ggplot(eigenvec_plot %>% filter(Species == "Pspp"), aes(x = PC3, y = PC5, color = Species)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  scale_shape_manual(values = c(0:10)) +
  xlab(paste0("PC3: ", round(percent_var[3], 2), "% variance")) +
  ylab(paste0("PC5: ", round(percent_var[5], 2), "% variance")) +
  geom_text_repel(data = subset(eigenvec_plot, Species == "Pspp"), 
                  aes(label = Geno_ID), size = 2, max.overlaps = Inf) +
  theme_bw()

# PCs 3 and 5 seem to do the best job of separating out the Pspp samples from the rest

# faceted by location
PC35_location <- ggplot(eigenvec_plot, aes(x = PC3, y = PC5, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC3: ", round(percent_var[3], 2), "% variance")) +
  ylab(paste0("PC5: ", round(percent_var[5], 2), "% variance")) +
  geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Location) +
  theme_bw()
PC35_location


# Check remaining PCs
# PC 6/7
ggplot(eigenvec_plot, aes(x = PC6, y = PC7, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC6: ", round(percent_var[6], 2), "% variance")) +
  ylab(paste0("PC7: ", round(percent_var[7], 2), "% variance")) +
  geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Location) +
  theme_bw()

# PC 8/9
ggplot(eigenvec_plot, aes(x = PC8, y = PC9, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC8: ", round(percent_var[8], 2), "% variance")) +
  ylab(paste0("PC9: ", round(percent_var[9], 2), "% variance")) +
  geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Location) +
  theme_bw()

# PC 9/10
ggplot(eigenvec_plot, aes(x = PC8, y = PC10, shape = Species, color = Location)) +
  scale_color_manual(values = mypal) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC9: ", round(percent_var[9], 2), "% variance")) +
  ylab(paste0("PC10: ", round(percent_var[10], 2), "% variance")) +
  geom_text_repel(aes(label = Genotype), size = 2, max.overlaps = Inf) +
  facet_wrap(~ Location) +
  theme_bw()
```
![image](https://github.com/user-attachments/assets/109cd80c-259e-4723-a5ce-b9c059116d41)
![image](https://github.com/user-attachments/assets/1a192a30-a429-4669-b439-7bf01d96750e)
![image](https://github.com/user-attachments/assets/991f224e-8a50-427a-83af-518325ec0770)



![image](https://github.com/user-attachments/assets/761a8c48-65be-4333-b37f-776d3d129a45)
![image](https://github.com/user-attachments/assets/72b7f514-8ab8-49a1-a130-4ab11b331a3c)
![image](https://github.com/user-attachments/assets/8aac0058-b0eb-4064-be11-85c8117c041d)
\  



Insights:

- PC 1 (15.8%) seems to separate out a cryptic species, present at Maloata and Vatia, that we did not visually distinguish in the field (we did NOT label them as Pspp).
    - 5 individuals from Vatia and many individuals from Fagasa fall at an intermediate position along this PC axis (potential hybrids?).
- PC 2 (12.6%) seems to separate out another cryptic species that we did not visually distinguish in the field, but which was only sampled at Ofu Pool 600.
    - No intermediate individuals between this group and the main cluster
    - To a lesser extent, this axis also separates out 6 of the 13 "Pspp" samples at it's opposite end, but many samples labeled as "Pver" also fall out at or near this position.
- PC 3 (11.2%) clearly separates out 6 of the 13 unknown "Pspp" samples we collected, in addition to a few samples (ALOF 01, MALO 21, and VATI 18) that were labeled as "Pver". Additionally, ALOF X1 falls at an intermediate position along this axis.
- PC 4 (10.4%) seems to separate out yet another cryptic species that we did not visually distinguish in the field, but which was only sampled at Faga'alu (13 individuals). Furthermore, it also splits the individuals from AOAA into two groups with 3 samples intermediate between them. In other words, this PC splits all individuals collected into 3 major clusters.
- PC 5 (9.8%) clusters individuals similarly to PC 3 (separates the same 6 "Pspp" samples, ALOF X1, and the three samples labeled as "Pver"), but also splits the other "Pver" samples from most locations into two groups (but the location of these groups in PC spaces varies by location). The exceptions to this were Fagatele and Leone, which had very little spread along this axis.


## Admixture/ancestry analysis
### Using ADMIXTURE v1.3.0

First need to convert PLINK2 files to PLINK1 .bed format:
```bash
crun.plink plink2 \
  --pfile pver_all_ld_pruned_0.2_genotypes \
  --make-bed \
  --out pver_all_ld_pruned_0.2_genotypes
```

### Run ADMIXTURE with cross-validation for K=1-10 across 3 replicate runs
First make lookup table for value of K and replicate number for each planned run:
```
1       1
1       2
1       3
2       1
2       2
2       3
3       1
3       2
3       3
4       1
4       2
4       3
5       1
5       2
5       3
6       1
6       2
6       3
7       1
7       2
7       3
8       1
8       2
8       3
9       1
9       2
9       3
10      1
10      2
10      3
```

Next, chromosome/scaffold names need to be changed in the .bim file becuase ADMIXTURE only allows them to be strings of numbers:
```bash
# create new bim file with altered chromosome names
sed -E 's/N[A-Z]_([0-9]+)\.1_Pverrucosa/\1/g' pver_all_ld_pruned_0.2_genotypes.bim > admixture.bim

# backup original bim and replace with new bim file
mv pver_all_ld_pruned_0.2_genotypes.bim pver_all_ld_pruned_0.2_genotypes_original.bim
mv admixture.bim pver_all_ld_pruned_0.2_genotypes.bim
```

`admixture_array.slurm`
```bash
#!/bin/bash
#SBATCH --job-name=admixture_array_pver_2025-07-14
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-30%30     # 10 K values * 3 replicates each = 30 tasks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=5-00:00:00

# Load modules
module load container_env admixture/1.3


# Set working directory
BASEDIR=/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/
cd $BASEDIR/admixture

# Define max K and replicates
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASEDIR/admixture/k_rep_table.txt)    # pulls out line from lookup table
read K REP <<< "$LINE"                                                         # uses the "here string" bash operator to provide a string as input to a command as if it came from stdin
                                                                               # "read" splits the string by whitespace and assigns each part to the variables "K" and "REP"

# Print info for logging
echo "Running ADMIXTURE for K=${K}, Replicate=${REP}, Task ID=${SLURM_ARRAY_TASK_ID}"

# Define input file
INPUT=$BASEDIR/vcf/pver_all_ld_pruned_0.2_genotypes.bed

# Output prefix (can help organize results)
OUTPUT_PREFIX="admixture_K${K}_rep${REP}"


# Run ADMIXTURE with 8 threads
crun.admixture admixture --seed=$SLURM_ARRAY_TASK_ID --cv=10 -j8 $INPUT $K > ${OUTPUT_PREFIX}.log 2>&1


# Rename output files to avoid overwriting
mv pver_all_ld_pruned_0.2_genotypes.${K}.Q ${OUTPUT_PREFIX}.Q
mv pver_all_ld_pruned_0.2_genotypes.${K}.P ${OUTPUT_PREFIX}.P

echo "Run $OUTPUT_PREFIX completed."
```

Compile cross-validation error results and log-likelihood across runs:
```bash
for FILE in admixture_K*_rep*.log; do
  K=$(echo $FILE | grep -oP 'K\d+')
  REP=$(echo $FILE | grep -oP 'rep\d+')
  CVERR=$(grep "CV error" $FILE | awk -F": " '{print $2}')
  LOGLIK=$(grep -oP 'Loglikelihood:\s+\K-?(.\d+\.?\d+?)$' "$FILE")
  echo -e "$K\t$REP\t$CVERR\t$LOGLIK"
done > cv_loglik_summary.tsv
```

Plot cross validation error to look for best K:
`plot_admixture_results.R`
```r
# Ancestry analysis with ADMIXTURE - Pver all samples
# 2025-07-14
# Jason A. Toy

rm(list = ls())

library(tidyverse)


setwd("/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/admixture")


# load cross validation errors for all runs (3 runs x 10 K values)
cv <- read_tsv("cv_summary.tsv", col_names = c("K", "Rep", "CV_Error")) %>%
  mutate(K = str_remove(K, "K") %>% as.numeric(),
         Rep = str_remove(Rep, "rep") %>% as.numeric)

# calculate mean and sd
cvsum <- cv %>% group_by(K) %>% summarize(Mean_CV = mean(CV_Error), SD = sd(CV_Error))

# plot CV Error vs. K
ggplot(cvsum, aes(x = K, y = Mean_CV)) +
  geom_line(color = "deepskyblue3", linewidth = 1) +
  geom_point(size = 2, color = "deepskyblue3") +
  geom_errorbar(aes(ymin = Mean_CV - SD, ymax = Mean_CV + SD), width = 0.1, color = "gray50") +
  scale_x_continuous(breaks = cvsum$K) +
  labs(title = "ADMIXTURE CV Error by K",
       x = "Number of Ancestral Populations (K)",
       y = "Mean CV Error ± SD") +
  theme_bw(base_size = 14)
```


There really isn't an "elbow" in the plot. The biggest decrease in CV is from K=1 to K=2 and the rest of the decreases are pretty equal and not that much smaller than the initial drop. One possiblility is that the best K is actually >10 so I created a new k_rep_table.txt file for K=11-20 and reran the `admixture_array.slurm` script.

Re-compile cross-validation error results across runs:
```bash
for FILE in admixture_K*_rep*.log; do
  K=$(echo $FILE | grep -oP 'K\d+')
  REP=$(echo $FILE | grep -oP 'rep\d+')
  CVERR=$(grep "CV error" $FILE | awk -F": " '{print $2}')
  echo -e "$K\t$REP\t$CVERR"
done > cv_summary.tsv
```

Re-plot cross validation error to look for best K:
`plot_admixture_results.R`
```r

```


## nMDS
Calculate Identity-by-state distance (1 - proportion of alleles shared) with PLINK1:
```bash
# load PLINK v1.9
module load plink/1.9-20240319

# calculate distance matrix
crun.plink plink \
	--bfile pver_all_ld_pruned_0.2_genotypes \
	--allow-extra-chr \
	--distance square 1-ibs \
	--out pver_all_ld_pruned_0.2_mdist


# calculate clustering and mds that may be useful for later analyses
crun.plink plink \
  --bfile pver_all_ld_pruned_0.2_genotypes \
  --distance square ibs \
  --cluster \
  --mds-plot 10 \
  --out pver_all_ld_pruned_0.2_ibsdist \
  --allow-extra-chr
```
parameters:

  --distance square ibs                       # calculates IBS distances
  --cluster                                   # does hierarchal clustering (UPGMA) for you based on the distance matrix
  --mds-plot 10                               # does MDS (PCoA) for you and gives coordinates as output
  --allow-extra-chr                           # allows nonstandard chromosome names


Run nMDS in R:
```r
library(vegan)

# Read distance matrix
dist_matrix <- as.matrix(read.table("pver_all_ld_pruned_0.2_mdist.mdist"))
sample_ids <- read.table("pver_all_ld_pruned_0.2_mdist.mdist.id")$V2

# Assign row/col names to distance matrix
rownames(dist_matrix) <- sample_ids
colnames(dist_matrix) <- sample_ids

# Compute nMDS
set.seed(123)  # for reproducibility
nmds_results <- metaMDS(dist_matrix, k = 2, trymax = 100)

# View stress
nmds_results$stress
  # stress is small (0.059; less than 0.1), which is good

# Get the scores (site coordinates in reduced space)
nmds_scores <- as.data.frame(scores(nmds_results))
nmds_scores$SampleID <- rownames(nmds_scores)

nmds_plotdat <- nmds_scores %>% 
  separate(SampleID, into = c("Year", "Location", "Species", "Genotype", "Rep"), sep = "_", remove = FALSE) %>% 
  mutate(Location = as.factor(Location),
         Species = as.factor(Species))

# Plot all samples
ggplot(nmds_plotdat, aes(x = NMDS1, y = NMDS2, shape = Species, color = Location)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = mypal) +
  scale_shape_manual(values = c(15, 16)) +
  geom_point(size = 2, alpha = 0.3) +
  theme_minimal() +
  labs(title = "nMDS all samples", x = "NMDS1", y = "NMDS2")

# Restrict x-axis to remove likely Pocillopora eydouxi
ggplot(nmds_plotdat, aes(x = NMDS1, y = NMDS2, shape = Species, color = Location)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = mypal) +
  scale_shape_manual(values = c(15, 16)) +
  geom_point(size = 2, alpha = 0.3) +
  xlim(c(0.0, 0.015)) +
  theme_minimal() +
  labs(title = "nMDS of Whole Genome SNP Data", x = "NMDS1", y = "NMDS2")
```

![image](https://github.com/user-attachments/assets/ccd71076-f579-4a20-911a-d42b9e07fbf1)
![image](https://github.com/user-attachments/assets/c93f432e-22ac-4539-8d30-6a752a26f286)


## Identify clones by recalculating 1-IBS distances based on unpruned SNP set
This is the set of 1,134,549 SNPs produced prior to LD-pruning, but after filtering by missingness, MAF, and allele number

First convert the VCF to PLINK binary format:
```bash
# Load PLINK2
module load plink/2024.03.02

crun.plink plink2 \
  --pfile pver_all_MISSMAFfiltered_genotypes \
  --make-bed \
  --out pver_all_MISSMAFfiltered_genotypes
```

Then calculate 1 - IBS distance matrix:
```bash
# Load PLINK1.9
module load plink/1.9-20240319

crun.plink plink\
  --bfile pver_all_MISSMAFfiltered_genotypes \
  --distance square 1-ibs \
  --out pver_all_MISSMAFfiltered_genotypes_mdist \
  --allow-extra-chr
```

Plot dendrogram in R:
```
# Read distance matrix (1-IBS) for UNPRUNED SNP set
dist_matrix_up <- as.matrix(read.table("pver_all_MISSMAFfiltered_genotypes_mdist.mdist"))
sample_ids_up <- read.table("pver_all_MISSMAFfiltered_genotypes_mdist.mdist.id")$V2

# Assign row/col names to distance matrix
rownames(dist_matrix_up) <- sample_ids_up
colnames(dist_matrix_up) <- sample_ids_up

# Convert distance matrix to a dist object
dist_obj_up <- as.dist(dist_matrix_up)

# Hierarchical clustering (default method = "complete")
hc_up <- hclust(dist_obj_up, method = "ward.D2")

# Simple plot
plot(hc_up, main = "Hierarchical Clustering of Samples (unpruned SNP set)", cex = 0.4)
```

Results aren't much different than the hierarchical clustering of the pruned SNP set. Let's try using identity-by-descent (IBD) instead.

## Identify clones by calculateing IBD with pruned dataset

Calculate proportion identity-by-descent (PI_HAT) with PLINK1 using --genome flag:
```bash
# load PLINK v1.9
module load plink/1.9-20240319

# calculate IBD (PI_HAT)
crun.plink plink \
	--bfile pver_all_ld_pruned_0.2_genotypes \
	--allow-extra-chr \
	--genome \
	--out pver_all_ld_pruned_0.2_relatedness
```

Now identify clones in R with a PI_HAT threshold:
```r
relatedness <- read.table("pver_all_ld_pruned_0.2_relatedness.genome", header = TRUE) %>% 
  separate(IID1, into = c("Year", "Location", "Species", "Genotype", "TechRep"), remove = FALSE) %>% 
  separate(IID2, into = c("Year2", "Location2", "Species2", "Genotype2", "TechRep2"), remove = FALSE) %>% 
  mutate(Sample1 = paste(Year, Location, Species, Genotype, sep = "_"),
         Sample2 = paste(Year2, Location2, Species2, Genotype2, sep = "_")
         )

techreps <- relatedness %>% filter(Sample1 == Sample2)
# Comparisons between the 4 samples with technical replicates have PI_HAT values of 0.8909, 0.8907, 0.9150, and 0.8998. This is lower than expected. Even with sequencing artifacts they should be > 0.95

# Plot distribution of PI_HAT values
ggplot(relatedness, aes(x = PI_HAT)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue") +
  theme_minimal() +
  labs(title = "Pairwise PI_HAT Distribution", x = "PI_HAT", y = "Count")
# Vast majority of comparisons have PI_HAT of 0

# Remove the 0 values to better visualize the rest of the distribution
ggplot(relatedness %>% filter(PI_HAT > 0), aes(x = PI_HAT)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue") +
  theme_minimal() +
  labs(title = "Pairwise PI_HAT Distribution (non-zero)", x = "PI_HAT", y = "Count")
```
Clones have PI_HAT values of 0.8909, 0.8907, 0.9150, and 0.8998. This is lower than expected. Even with sequencing artifacts they should be > 0.95


Compute genotype discordance with bcftools:
```bash
crun.bcftools bcftools gtcheck -p 2024_LEON_Pver_14_1,2024_LEON_Pver_14_2 pver_all_ld_pruned_0.2_genotypes.vcf > 2024_LEON_Pver_14.gtcheck
crun.bcftools bcftools gtcheck -p 2024_OFU6_Pver_09_1,2024_OFU6_Pver_09_2 pver_all_ld_pruned_0.2_genotypes.vcf > 2024_OFU6_Pver_09.gtcheck
crun.bcftools bcftools gtcheck -p 2024_FALU_Pver_23_1,2024_FALU_Pver_23_2 pver_all_ld_pruned_0.2_genotypes.vcf > 2024_FALU_Pver_23.gtcheck
crun.bcftools bcftools gtcheck -p 2024_FASA_Pver_34_1,2024_FASA_Pver_34_2 pver_all_ld_pruned_0.2_genotypes.vcf > 2024_FASA_Pver_34.gtcheck

# add header line and combine across samples
{
grep "^#DC" 2024_LEON_Pver_14.gtcheck;
for FILE in *.gtcheck; do
	tail -n 1 $FILE
done;
} > all_techreps.gtchecksummary
```
```
#DC     [2]Query Sample [3]Genotyped Sample     [4]Discordance  [5]-log P(HWE)  [6]Number of sites compared
DC      2024_FALU_Pver_23_1     2024_FALU_Pver_23_2     1.790155e+04    2.404484e+05    80161
DC      2024_FASA_Pver_34_1     2024_FASA_Pver_34_2     1.826797e+04    2.421047e+05    80154
DC      2024_LEON_Pver_14_1     2024_LEON_Pver_14_2     1.575252e+04    2.428629e+05    80172
DC      2024_OFU6_Pver_09_1     2024_OFU6_Pver_09_2     1.707957e+04    2.411082e+05    80173
```
These correspond to mismatch rates of ~20%, which is a lot more than expected for technical replicates.


To see if SNP filtering is creating this issue rerun with less filtered VCFs:
```bash
crun.bcftools bcftools gtcheck -p 2024_LEON_Pver_14_1,2024_LEON_Pver_14_2 pver_all_MISSMAFfiltered_genotypes.vcf > 2024_LEON_Pver_14_unpruned.gtcheck
crun.bcftools bcftools gtcheck -p 2024_OFU6_Pver_09_1,2024_OFU6_Pver_09_2 pver_all_MISSMAFfiltered_genotypes.vcf > 2024_OFU6_Pver_09_unpruned.gtcheck
crun.bcftools bcftools gtcheck -p 2024_FALU_Pver_23_1,2024_FALU_Pver_23_2 pver_all_MISSMAFfiltered_genotypes.vcf > 2024_FALU_Pver_23_unpruned.gtcheck
crun.bcftools bcftools gtcheck -p 2024_FASA_Pver_34_1,2024_FASA_Pver_34_2 pver_all_MISSMAFfiltered_genotypes.vcf > 2024_FASA_Pver_34_unpruned.gtcheck

# add header line and combine across samples
{
grep "^#DC" 2024_LEON_Pver_14_unpruned.gtcheck;
for FILE in *_unpruned.gtcheck; do
	tail -n 1 $FILE
done;
} > all_techreps_unpruned.gtchecksummary
```
```
#DC     [2]Query Sample [3]Genotyped Sample     [4]Discordance  [5]-log P(HWE)  [6]Number of sites compared
DC      2024_FALU_Pver_23_1     2024_FALU_Pver_23_2     4.974384e+04    3.459043e+06    1134466
DC      2024_FASA_Pver_34_1     2024_FASA_Pver_34_2     5.648805e+04    3.453628e+06    1134454
DC      2024_LEON_Pver_14_1     2024_LEON_Pver_14_2     4.232621e+04    3.481962e+06    1134487
DC      2024_OFU6_Pver_09_1     2024_OFU6_Pver_09_2     4.966462e+04    3.452874e+06    1134489
```
This results in much lower discordance rates of ~5%.


Recalculate again using the quality-filtered only vcf (5.255 million SNPs) to get rid of potential issues with MAF-filtering, which may also affect PI_HAT calculations:
```bash
crun.bcftools bcftools gtcheck -p 2024_LEON_Pver_14_1,2024_LEON_Pver_14_2 pver_all_QDPfiltered_genotypes.vcf.gz > 2024_LEON_Pver_14_noPLINKfilt.gtcheck
crun.bcftools bcftools gtcheck -p 2024_OFU6_Pver_09_1,2024_OFU6_Pver_09_2 pver_all_QDPfiltered_genotypes.vcf.gz > 2024_OFU6_Pver_09_noPLINKfilt.gtcheck
crun.bcftools bcftools gtcheck -p 2024_FALU_Pver_23_1,2024_FALU_Pver_23_2 pver_all_QDPfiltered_genotypes.vcf.gz > 2024_FALU_Pver_23_noPLINKfilt.gtcheck
crun.bcftools bcftools gtcheck -p 2024_FASA_Pver_34_1,2024_FASA_Pver_34_2 pver_all_QDPfiltered_genotypes.vcf.gz > 2024_FASA_Pver_34_noPLINKfilt.gtcheck

# add header line and combine across samples
{
grep "^#DC" 2024_LEON_Pver_14_noPLINKfilt.gtcheck;
for FILE in *_noPLINKfilt.gtcheck; do
	tail -n 1 $FILE
done;
} > all_techreps_noPLINKfilt.gtchecksummary
```
```
#DC     [2]Query Sample [3]Genotyped Sample     [4]Discordance  [5]-log P(HWE)  [6]Number of sites compared
DC      2024_FALU_Pver_23_1     2024_FALU_Pver_23_2     3.751238e+05    4.435625e+07    4850604
DC      2024_FASA_Pver_34_1     2024_FASA_Pver_34_2     3.918278e+05    4.433923e+07    4850604
DC      2024_LEON_Pver_14_1     2024_LEON_Pver_14_2     3.689223e+05    4.440351e+07    4850604
DC      2024_OFU6_Pver_09_1     2024_OFU6_Pver_09_2     3.564821e+05    4.442604e+07    4850604
```
This results in discordance rates of 7-8%.



MAF filtering can also decrease PI_HAT estimates between technical replicates because it removes additional sites where replicates would be identical.
It also enriches for common variants which are more likely to be heterozygous. Heterozygous sites are more prone to genotyping error or allelic dropout, which can also deflate PI_HAT between replicates.

So lets try refiltering the initial quality/depth-filtered VCF without a MAF filter (filter based on missingess and remove indels and multiallelic SNPs):
```bash
crun.plink plink2 \
  --vcf pver_all_QDPfiltered_genotypes.vcf.gz \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno 0.2 \
  --mind 0.2 \
  --make-pgen \
  --out pver_all_MISSfiltered_genotypes
```
```
...
0 samples removed due to missing genotype data (--mind).
396 samples (0 females, 0 males, 396 ambiguous; 396 founders) remaining after main filters.
Calculating allele frequencies... done.
--geno: 8 variants removed due to missing genotype data.
4585661 variants remaining after main filters.
```


Or with a much lower MAF filter:
```bash
module load plink/2024.03.02    # Note: this is PLINK2

crun.plink plink2 \
  --vcf pver_all_QDPfiltered_genotypes.vcf.gz \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno 0.2 \
  --mind 0.2 \
  --maf 0.001 \
  --make-pgen \
  --out pver_all_MISSMAF001filtered_genotypes


  # geno 0.2    removes SNPs where more than 20% of individuals have missing genotypes
  # mind 0.2    removes individuals who are missing more than 20% of genotype data across all SNPs
  # maf 0.05    removes SNPs where the minor allele frequency is less than 0.05
  # threads N    multithreading option available
```
```
...
0 samples removed due to missing genotype data (--mind).
396 samples (0 females, 0 males, 396 ambiguous; 396 founders) remaining after main filters.
Calculating allele frequencies... done.
--geno: 8 variants removed due to missing genotype data.
42988 variants removed due to allele frequency threshold(s)
4542673 variants remaining after main filters.
```
The MAF 0.001 filter only removes 42,988 SNPs, so I will continue with this set for recalculating PI_HAT

Convert to VCF for viewing:
```bash
crun.plink plink2 \
  --pgen pver_all_MISSMAF001filtered_genotypes.pgen \
  --psam pver_all_MISSMAF001filtered_genotypes.psam \
  --pvar pver_all_MISSMAF001filtered_genotypes.pvar \
  --recode vcf \
  --out pver_all_MISSMAF001filtered_genotypes
```
This vcf has **4,542,673** SNPs remaining.

Recalculate discordance between technical replicates
```bash
crun.bcftools bcftools gtcheck -p 2024_LEON_Pver_14_1,2024_LEON_Pver_14_2 pver_all_MISSMAF001filtered_genotypes.vcf > 2024_LEON_Pver_14_MAF001unpruned.gtcheck
crun.bcftools bcftools gtcheck -p 2024_OFU6_Pver_09_1,2024_OFU6_Pver_09_2 pver_all_MISSMAF001filtered_genotypes.vcf > 2024_OFU6_Pver_09_MAF001unpruned.gtcheck
crun.bcftools bcftools gtcheck -p 2024_FALU_Pver_23_1,2024_FALU_Pver_23_2 pver_all_MISSMAF001filtered_genotypes.vcf > 2024_FALU_Pver_23_MAF001unpruned.gtcheck
crun.bcftools bcftools gtcheck -p 2024_FASA_Pver_34_1,2024_FASA_Pver_34_2 pver_all_MISSMAF001filtered_genotypes.vcf > 2024_FASA_Pver_34_MAF001unpruned.gtcheck

# add header line and combine across samples
{
grep "^#DC" 2024_LEON_Pver_14_MAF001unpruned.gtcheck;
for FILE in *_MAF001unpruned.gtcheck; do
	tail -n 1 $FILE
done;
} > all_techreps_MAF001unpruned.gtchecksummary
```
```
#DC     [2]Query Sample [3]Genotyped Sample     [4]Discordance  [5]-log P(HWE)  [6]Number of sites compared
DC      2024_FALU_Pver_23_1     2024_FALU_Pver_23_2     7.220304e+04    4.192155e+07    4542558
DC      2024_FASA_Pver_34_1     2024_FASA_Pver_34_2     7.927406e+04    4.193757e+07    4542543
DC      2024_LEON_Pver_14_1     2024_LEON_Pver_14_2     6.182430e+04    4.196157e+07    4542588
DC      2024_OFU6_Pver_09_1     2024_OFU6_Pver_09_2     7.100474e+04    4.192378e+07    4542581
```
This results in discordance rates between 1.3% and 1.8%. This is much more reasonable for technical replicates.


To see how this compares to the similarity between other samples, recalculate PI_HAT:
Convert VCF to PLINK1 binary format:
```bash
# Load PLINK2
module load plink/2024.03.02

crun.plink plink2 \
  --pfile pver_all_MISSMAF001filtered_genotypes \
  --make-bed \
  --out pver_all_MISSMAF001filtered_genotypes
```

Run calculation:
```bash
# load PLINK v1.9
module load plink/1.9-20240319

# calculate IBD (PI_HAT)
crun.plink plink \
	--bfile pver_all_MISSMAF001filtered_genotypes \
	--allow-extra-chr \
	--genome full \
	--out pver_all_MISSMAF001filtered_genotypes_relatedness
```
Relatedness (PI_HAT) of technical replicates is higher, but so is the relatedness of all other sample comparisons, so we still have the same issue.


PI_HAT assumes an outbred population, but Pocillopora spp are known for high clonality and brooding of larvae, so this assumption does not hold. Instead try calculating relatedness with VCFtools.
This method doesn't assume distributions based on HWE, but instead uses allele frequencies from your data, so is less sensitive to inbreeding.
```bash
#First downgrade VCF from v4.3 to v4.2 so it will work with VCFtools
sed 's/VCFv4.3/VCFv4.2/' pver_all_MISSMAF001filtered_genotypes.vcf | crun.bcftools bgzip > pver_all_MISSMAF001filtered_genotypes_v42.vcf.gz


module load dosage_convertor
crun.dosage_convertor vcftools --gzvcf pver_all_MISSMAF001filtered_genotypes_v42.vcf.gz --relatedness2 --out pver_all_MISSMAF001filtered_genotypes_v42.relatedness2
```

Visualize results in R













Take a closer look at quality and depth of discordant sites for a given sample comparison.

First compress vcfs into vcf.gz with bgzip:
```bash
# Compress the VCF file with bgzip
crun.bcftools bgzip -c pver_all_MISSMAFfiltered_genotypes.vcf > pver_all_MISSMAFfiltered_genotypes.vcf.gz
crun.bcftools bgzip -c pver_all_ld_pruned_0.2_genotypes.vcf > pver_all_ld_pruned_0.2_genotypes.vcf.gz

# Index the compressed file
crun.bcftools bcftools index pver_all_MISSMAFfiltered_genotypes.vcf.gz
crun.bcftools bcftools index pver_all_ld_pruned_0.2_genotypes.vcf.gz
```


```bash
#!/bin/bash

# Set variables
VCF="pver_all_MISSMAFfiltered_genotypes.vcf"
SAMPLE1="2024_LEON_Pver_14_1"
SAMPLE2="2024_LEON_Pver_14_2"
OUT="discordant_sites_${SAMPLE1}_vs_${SAMPLE2}.tsv"

# Extract discordant sites
crun.bcftools bcftools gtcheck "$VCF" -p ${SAMPLE1},${SAMPLE2} -v discordant | \
    awk '!/^#/ {print $1"\t"$2}' > discordant_positions.txt

# Create region list for bcftools view
awk '{print $1":"$2"-"$2}' discordant_positions.txt > discordant_regions.txt

# Extract GT, DP, GQ for both samples at discordant sites
crun.bcftools bcftools view -R discordant_regions.txt -s ${SAMPLE1},${SAMPLE2} -f PASS "$VCF" | \
    crun.bcftools bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t%DP\t%GQ\t]\n' > "$OUT"

echo "Done. Output written to $OUT"
```

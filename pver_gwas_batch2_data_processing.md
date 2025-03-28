Processing of 2024 P. verrucosa WGS data from American Samoa using GATK4 tools and adapted best practices
Data Batch 2
================
Jason Toy
2025-03-27

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


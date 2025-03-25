Pver GWAS Pilot Data Processing
================
Jason Toy
2024-07-18

# Part I - Data Processing

## Setup sample table and sample list

 

Based on the Therkildsen Lab pipeline, this pipeline is set up to link
up fastq prefix names with sample details based on a **fastq table** in
`$BASEDIR/sample_lists/fastq_table_pver_pilot.txt`, set up as described
below. Create these files for your samples before starting your
analysis.

For the scripts below to work, the sample table has to be a tab
deliminated table with the following six columns, strictly in this
order:

- `prefix` the prefix of raw fastq file names

- `lane_number` lane number; each sequencing lane or batch should be
  assigned a unique identifier. This is important so that if you
  sequence a library across multiple different sequencing lanes, you can
  keep track of which lane/batch a particular set of reads came from
  (important for accounting for sequencing error patterns or batch
  effects).

- `seq_id` sequence ID; this variable is only relevant when different
  libraries were prepared out of the same sample and were run in the
  same lane (e.g. if you wanted to include a replicate). In this case,
  seq_id should be used to distinguish these separate libraries. If you
  only have a single library prepared from each of your samples (even if
  you sequence that same library across multiple lanes), you can just
  put 1 for all samples in this column.

- `sample_id` sample ID; a unique identifier for each individual
  sequenced

- `population` population name; the population or other relevant
  grouping variable that the individual belongs to

- `data_type` data type; there are only two allowed entries here: pe
  (for paired-end data) or se (for single end data). We need this in the
  table because for some of our processing steps, the commands are
  slightly different for paired-end and single-end data.

It is important to make sure that the combination of sample_id, seq_id,
and lane_number is unique for each fastq file.

You can also add additional columns to the end of the table as you see
fit, but the first 6 must be in the order listed above.  

We’ll also use a second file that we call a **fastq list**. This is
simply a list of prefixes for the samples we want to analyze. Our sample
table can contain data for all individuals in our study, but at any
given time, we may only want to perform an operation on a subset of
them. Note that it’s just a list of fastq name prefixes, each on a
separate line and there should be no header in this file. If the samples
to be analyzed include all the samples in the **fastq table**, the
**fastq list** can be created by running the following command:

``` bash
cut -f1 fastq_table_pver_pilot.txt | tail -n +2 > fastq_list_pver_pilot.txt
```

 

## Assign starting environment variables

### Then check that the environment is working as expected by using a for loop to print the sample IDs for all sample files in your sample list:

``` bash
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/
# Change this to your base working directory

SAMPLELIST=$BASEDIR/sample_lists/fastq_list_pver_pilot.txt
# Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the fastq table.

SAMPLETABLE=$BASEDIR/sample_lists/fastq_table_pver_pilot.txt
# Path to a fastq table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID.




#Use a for loop to print the sample IDs for all sample files in your sample list:

for SAMPLEFILE in `cat $SAMPLELIST`; do   # Loop through each of the prefixes listed in our fastq list

    # For each prefix, extract the associated sample ID (column 4) and population (column 5) from the table
    SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
    POPULATION=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
    echo $SAMPLEFILE refers to sample $SAMPLE_ID from $POPULATION

done
```

 

## Check read quality and trim reads using fastp (v0.23.2)

Make a new directory in the \$BASEDIR/reports directory called
`fastp_reports`

``` bash
mkdir $BASEDIR/reports/fastp_reports
```

 

Create array SLURM script to run fastp on all samples at once:

``` bash
nano fastp_array.slurm
```

``` bash
#!/bin/bash

#SBATCH --job-name fastp_array_2024-07-22
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env
module load fastp
                                                                                                                                     ## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy/
RAWDATA=$BASEDIR/raw_sequence_data/2024-07-11_TxGenPverrucosaTestSet/24200Brs_N24114 #path to raw fq.gz files
OUTDIR=$BASEDIR/pver_gwas_pilot/trimmed_fastq
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_list_pver_pilot.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_table_pver_pilot.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
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
SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely

echo $SAMPLE_UNIQ_ID

## Define the output path and file prefix
SAMPLETRIM=$OUTDIR/$SAMPLE_UNIQ_ID


## Run fastp
crun.fastp fastp -i $RAWDATA/$SAMPLEFILE$RAW_R1 \
    -I $RAWDATA/$SAMPLEFILE$RAW_R2 \
    -o ${SAMPLETRIM}_f_paired_trim.fastq.gz \
    -O ${SAMPLETRIM}_r_paired_trim.fastq.gz \
    --adapter_fasta $BASEDIR/pver_gwas_pilot/adapters.fa \
    --cut_tail \
    -l 40 \
    -h $BASEDIR/pver_gwas_pilot/reports/fastp_reports/${SAMPLE_UNIQ_ID}_fastp.html \
    -j $BASEDIR/pver_gwas_pilot/reports/fastp_reports/${SAMPLE_UNIQ_ID}_fastp.json \
    --thread 20
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

crun.multiqc multiqc --interactive --filename multiqc_report_pver_pilot_fastp_interactive .
# --interactive forces the creation of an interactive html file even when sample sizes are high (instead of a flat file)
```

 

### Run FastQC on trimmed (fastp) reads. First make new directory for FastQC reports and sample list for trimmed fastq files:

``` bash
# make new directory
mkdir $BASEDIR/reports/fastqc_trimmed_reports

# make new sample list
cd $BASEDIR/trimmed_fastq

ls *.fastq.gz > $BASEDIR/sample_lists/fastq_trimmed_list_pver_pilot.txt


# if you just want a list of single prefixes for each samples (instead of two file names (f and r) per sample), you can make this with the following commands:
for FILE in `ls *.fastq.gz`; do
  echo ${FILE%_r_*} >> $BASEDIR/sample_lists/temp.txt
done

grep -v "paired_trim.fastq.gz" $BASEDIR/sample_lists/temp.txt  > $BASEDIR/sample_lists/fastq_trimmed_list_pver_pilot_single_prefixes.txt

rm $BASEDIR/sample_lists/temp.txt
```

 

Then run array script `fastqc_array.slurm`

``` bash
#!/bin/bash

#SBATCH --job-name fastqc_array_2024-07-23
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-32%32
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env
module load fastqc/0.11.9

## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/
SEQDATA=$BASEDIR/trimmed_fastq/ #path to fq.gz files
OUTDIR=$BASEDIR/reports/fastqc_trimmed_reports/
SAMPLELIST=$BASEDIR/sample_lists/fastq_trimmed_list_pver_pilot.txt #Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Run FastQC
crun.fastqc fastqc -t 16 -o $OUTDIR $SEQDATA/$SAMPLEFILE
```

 

``` bash
sbatch $BASEDIR/scripts/fastqc_array.slurm
```

 

Compile with MultiQC:

``` bash
cd $BASEDIR/reports/fastqc_trimmed_reports/

module load container_env
module load multiqc/1.13

crun.multiqc multiqc --interactive --filename multiqc_report_pver_pilot_fastqc_trimmed_interactive .
# --interactive forces the creation of an interactive html file even when sample sizes are high (instead of a flat file)
```

 

### Trimming results for all samples:

The level of duplication estimated by fastp is quite different from
FastQC. This seems to be a known/expected result because FastQC does not
analyze read data in pairs, while fastp does. FastQC has the duplicate
percentage at around 40% both before and after trimming, while fastp
puts it closer to 20-25%. Even if we assume the true duplication level
is something like 20-25%, this is still a bit disappointing and
indicates a bit of a library issue. It seems like there are two possible
causes. One is that the libraries made by TxGen tend to be
low-complexity, or possibly they didn’t load enough library into the
machine after they decided to double our output and ended up with lots
of optical duplicates? Likely the former. This also is probably caused
by the amount of data we got. Any library will at some point be
sequenced to its maximum potential and additional sequencing will only
yield duplicates. Hopefully this will be less of an issue when we
multiplex many more samples on the same flowcell, but hard to say.

For reference, the fastp-estimated duplication level for the Ahya pilot
set we sent to Berkeley was 2-4%. My guess is that this is in part
because we ordered much much less sequencing from them for these
samples, but also probably because their library prep yields higher
complexity. The Ahya data we got from TxGen was closer to the 15-20%
range. I was hoping for better rates this time because our DNA was much
better quality.

### Subsetting Test

We can test the idea that higher duplicates are a result of a greater
amount of sequencing by subsetting the data and rerunning FastQC. The
`seqtk_subsetting_array.slurm` script takes 40 million random reads from
each file (80 million total per sample).

``` bash
#!/bin/bash

#SBATCH --job-name seqtk_subsetting_array_2024-07-24
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-32%32
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env
module load seqtk

## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/
SEQDATA=$BASEDIR/trimmed_fastq/ #path to fq.gz files
OUTDIR=$BASEDIR/trimmed_fastq/jt_subsetting
SAMPLELIST=$BASEDIR/sample_lists/fastq_trimmed_list_pver_pilot.txt #Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Run FastQC
crun.seqtk seqtk sample -s100 $SEQDATA/$SAMPLEFILE 40000000 | gzip > $OUTDIR/${SAMPLEFILE%.*}'subset40M.fastq.gz'

## Run FastQC
module load fastqc/0.11.9

crun.fastqc fastqc -t 16 -o $OUTDIR $OUTDIR/${SAMPLEFILE%.*}'subset40M.fastq.gz'
```

 

Move fastqc output files to new directory

``` bash
mkdir fastqc_reports_subset
mv *fastqc.* fastqc_reports_subset/ 
```

 

Run MultiQC on fastqc files:

``` bash
cd fastqc_reports_subset

module load container_env multiqc/1.13

crun.multiqc multiqc --interactive --filename multiqc_report_pver_pilot_fastqc_trimmed_subset80M_interactive .
```

 

### Results of subsetting test

Subsetting indeed resulted in a ~50% reduction in the percentage of
duplicates down to the 12-20% range (as caculated by FastQC). This
indicates to me that the libraries are not of high enough complexity to
be sequenced to this degree.  

I repeated this subsetting process at 270M reads per sample (135M per
file), and 150M reads per sample (75M per file) to see how duplication
rate decreases with decreased sequencing effort. Note that the data are
plotted for each read (F & R) separately:
![](pver_gwas_pilot_data_processing_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## ITS2 symbiont typing

Before we can map the reads, we need to know what symbiont clades are
present so that we can create a combined holobiont genome.

This time, I am using a version of the SymPortal database that Dan and I
edited to trim the end of the “I” clade sequences that were causing
erroneous matches. This file was uploaded to my `references` directory
and is called `Voolstra_SymPortal_published_div_20230315_Itrimmed.fasta`

Now we need to index the database so we can map to it:

``` bash
salloc --partition=main --exclusive
tmux

module load container_env
module load bowtie2/2.4.1

cd /cm/shared/courses/dbarshis/barshislab/jtoy/references

crun.bowtie2 bowtie2-build Voolstra_SymPortal_published_div_20230315_Itrimmed.fasta Voolstra_SymPortal_published_div_20230315_Itrimmed
#second argument is base name for index files that will be created
```

 

Then map the trimmed fastq files to the reference using the following
array job script, `its2_mapping_array.slurm`:

``` bash
#!/bin/bash

#SBATCH --job-name its2_mapping_array_2024-07-24
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env
module load bowtie2/2.4.1

## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
FASTQDIR=$BASEDIR/pver_gwas_pilot/trimmed_fastq #path to trimmed fastq.gz files
OUTDIR=$BASEDIR/pver_gwas_pilot/its2_mapping
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_list_pver_pilot.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_table_pver_pilot.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
FASTQ_SUFFIX_1=_f_paired_trim.fastq.gz # Suffix to trimmed fastq files. Use forward reads with paired-end data.
FASTQ_SUFFIX_2=_r_paired_trim.fastq.gz # Suffix to trimmed fastq files. Use reverse reads with paired-end data.
REFBASENAME=Voolstra_SymPortal_published_div_20230315_Itrimmed

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
SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2` # Define platform unit (PU), which is the lane number

echo $SAMPLE_UNIQ_ID

## Define the output path and file prefix
SAMPLEOUT=$OUTDIR/$SAMPLE_UNIQ_ID


## Run bowtie2
crun.bowtie2 bowtie2 -p 30 --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA \
        --local -x $BASEDIR/references/$REFBASENAME -1 $FASTQDIR/$SAMPLE_UNIQ_ID$FASTQ_SUFFIX_1 -2 $FASTQDIR/$SAMPLE_UNIQ_ID$FASTQ_SUFFIX_2 \
        -k 5 -S $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'


## Keep only best alignment (remove secondary alignments)
module unload bowtie2
module load container_env samtools
crun.samtools samtools view -@38 -h $SAMPLEOUT'_bt2_'$REFBASENAME'.sam' | grep -v "YT:Z:UP" | grep -vP "\s255\s" > $SAMPLEOUT'_bt2_'$REFBASENAME'_k1.sam'


# When the -k argument is used in bowtie2, the supplemental alignments (i.e. not the best alignment) are given a MAPQ score of 255. bowtie2 also reports unaligned reads in this mode with the YT:Z:UP code. To remove these unnecessary lines from the alignment file, the above samtools/grep command is used.
```

``` bash
sbatch its2_mapping_array.slurm
```

 

Query-sort, validate, then dedup alignments:

``` bash
#!/bin/bash

#SBATCH --job-name its2_dedup_array_2024-09-19
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env
module load bowtie2/2.4.1

## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
FASTQDIR=$BASEDIR/pver_gwas_pilot/trimmed_fastq #path to trimmed fastq.gz files
OUTDIR=$BASEDIR/pver_gwas_pilot/its2_mapping
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_list_pver_pilot.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_table_pver_pilot.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
FASTQ_SUFFIX_1=_f_paired_trim.fastq.gz # Suffix to trimmed fastq files. Use forward reads with paired-end data.
FASTQ_SUFFIX_2=_r_paired_trim.fastq.gz # Suffix to trimmed fastq files. Use reverse reads with paired-end data.
REFBASENAME=Voolstra_SymPortal_published_div_20230315_Itrimmed

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
SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2` # Define platform unit (PU), which is the lane number

echo $SAMPLE_UNIQ_ID

## Define the output path and file prefix
SAMPLEOUT=$OUTDIR/$SAMPLE_UNIQ_ID

#add version number to bowtie2 command header line in SAM file
sed 's/PN:bowtie2\tVN:/PN:bowtie2\tVN:2.4.1/' $SAMPLEOUT'_bt2_'$REFBASENAME'_k1.sam' > $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered.sam'

#remove SAM file
#rm $SAMPLEOUT'_bt2_'$REFBASENAME'_k1.sam'

## Change modules
module unload bowtie2
module load container_env gatk
GATK='crun.gatk gatk'

## Query-sort for duplicate removal with GATK
# Run SortSam to sort by query name and convert to BAM
$GATK --java-options "-Xmx100G" SortSam \
  --INPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered.sam' \
  --OUTPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted.bam' \
  --SORT_ORDER queryname

#sort and convert to BAM at the same time
#crun.samtools samtools sort -@ 40 -o $SAMPLEOUT'_bt2_'$REFBASENAME'_sorted.bam' $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'

# Run validation of BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted.val' \
  -M VERBOSE

#remove SAM file
rm $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered.sam'

## Mark and remove duplicates
$GATK --java-options "-Xmx100G" MarkDuplicates \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted_dedup.bam' \
  --METRICS_FILE $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted_dupstat.txt' \
  --REMOVE_DUPLICATES true

## Run validation of deduped BAM file
$GATK --java-options "-Xmx100G" ValidateSamFile \
  -I $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted_dedup.bam' \
  -O $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted_dedup.val' \
  -M VERBOSE

## Run SortSam to sort by coordinate for downstream processing
$GATK --java-options "-Xmx100G" SortSam \
  --INPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted_dedup.bam' \
  --OUTPUT $SAMPLEOUT'_bt2_'$REFBASENAME'_k1_reheadered_qsorted_dedup_coordsorted.bam' \
  --SORT_ORDER coordinate
  --CREATE_INDEX true

#index coordinate-sorted BAM with samtools
#module load container_env samtools
#crun.samtools samtools index -@ 40 $SAMPLEOUT'_bt2_'$REFBASENAME'_qsorted_dedup_coordsorted.bam'
```

Count reads mapping to each contig (strain) using
`countxpression_SB_advbioinf.py`. This will create a counts file for
every .sam file:

``` bash
module load container_env
module load python2

crun.python2 python2 $BASEDIR/pver_gwas_pilot/scripts/countxpression_SB_advbioinf.py $BASEDIR/pver_gwas_pilot/its2_mapping/*.sam
```

 

Compile counts files into one data file using
`ParseExpression2BigTable_advbioinf.py`.

First create list of contigs (strain names):

``` bash

cd /cm/shared/courses/dbarshis/barshislab/jtoy/references

grep ">" Voolstra_SymPortal_published_div_20230315_Itrimmed.fasta | cut -c 2- > its2_strain_names.txt

# Add header "sym" to list file
nano its2_strain_names.txt

cat its2_strain_names.txt
```

    sym
    C42.2
    C1
    C42i
    C1b
    C115k
    C3mk
    C42m
    C42bw
    C1j

 

Now run the python script

``` bash
cd $BASEDIR/its2_mapping

crun.python2 python2 ../scripts/ParseExpression2BigTable_advbioinf.py ../../references/its2_strain_names.txt pver_pilot_its2_counts_merged.txt no_match *20230315_Itrimmed_counts.txt

    # Script usage:
      #sys.argv[1] Input file name #list of your "gene" names (list must have a header line!)
      #sys.argv[2] Output file name
      #sys.argv[3] Text to add when no match
      #sys.argv[4:] Any number of files to add columns from
```

 

The output file shows the number of reads mapping to each strain with
one column for each sample:

    A1      0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
    A1.2    0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
    A10     0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
    A10a    0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
    A10b    0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
    A11     0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
    A12     0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
    A13     0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
    A13a    0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0

 

Start an interactive RStudio Server session or download results file to
local computer to run in Rstudio. Calculate proportion of reads for each
symbiont clade for each sample and plot them.

``` r
library(tidyverse)
library(reshape2)
library(RColorBrewer)


setwd("C:/Users/jason/OneDrive/Documents/odu/projects/global_search_genomics/pver_gwas/")

# load count data from ParseExpression2BigTable_advbioinf.py
raw <- read_tsv("pver_pilot_its2_counts_merged.txt", col_names = TRUE)

# make column names shorter
col_names <- c("sym", "VATI_21", "VATI_22", "VATI_23", "VATI_24", "VATI_25", "VATI_26", "VATI_27", "VATI_28", "VATI_29", "VATI_30", "VATI_31", "VATI_32", "VATI_33", "VATI_34", "VATI_35", "VATI_36")

colnames(raw) <- col_names

# reshape data to long format for summarizing and plotting
counts <- melt(raw) %>% mutate(clade = substr(sym,1,1)) %>% rename(sample = variable, count = value)

# summarize counts by sample and clade
sum_table <- counts %>% group_by(sample, clade) %>% summarise(clade_count = sum(count))

# plot total counts
ggplot(sum_table, aes(fill=clade, y=clade_count, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Set1") +
  ylab("read count") +
  ggtitle("Raw Read Counts - ITS2 Symbiont Mapping - P. verrucosa - Vatia 2024") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_blank())

ggsave(plot = last_plot(), file = "2024_pver_pilot_symbiont_total_counts_plot.png", width = 8, height = 8, units = "in")

# calculate percent reads for each clade for each sample
total_reads <- counts %>% group_by(sample) %>% summarise(total_count = sum(count))

# merge data frames to add percent column to summary table
sum_table_perc <- left_join(sum_table, total_reads, by = "sample") %>% mutate(perc = clade_count/total_count)

# plot read counts as percentages
ggplot(sum_table_perc, aes(fill=clade, y=perc, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Set1") +
  ylab("% of mapped reads") +
  ggtitle("Read Proportions - ITS2 Symbiont Mapping - P. verrucosa - Vatia 2024") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_blank())


ggsave(plot = last_plot(), file = "2024_pver_pilot_symbiont_proportions_plot.png", width = 8, height = 8, units = "in")
```

 

![](C:/Users/jason/OneDrive/Documents/odu/projects/global_search_genomics/pver_gwas/2024_pver_pilot_symbiont_total_counts_plot_old.png)

![](C:/Users/jason/OneDrive/Documents/odu/projects/global_search_genomics/pver_gwas/2024_pver_pilot_symbiont_proportions_plot_old.png)  

Results Summary: - Low overall mapping rates indicates that most of our
sequence data is from the host (which is normally a good thing) - Mostly
Durisdinium reads with some Cladocopium and a sprinkling of other clades
 

Convert .sam files to .bam files to save space:

``` bash
module load samtools

for FILE in *.sam; do
    BASENAME=`basename $FILE .sam`
    echo $BASENAME
    crun.samtools samtools view -Sb -@ 38 -O BAM -o $BASENAME.bam $FILE
done
```

 

### Compare counts from countexpression.py script to output summary from bowtie2

Run the following script (`parse_bowtie2_output_its2.sh`) to parse
summary mapping data from each file. The output file is called
`its2_bowtie2_mapping_summary.tsv`

``` bash
#!/bin/bash

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/

cd $BASEDIR/scripts/errfiles

ls *its2_mapping_array*.err > its2_errfile_list.txt

# create data file and add header if it is empty
header="sample\tconcordantly_0_times\tconcordantly_1_time\tconcordantly_2_or_more_times\toverall_rate" #assign header value
outfile="its2_bowtie2_mapping_summary.tsv"

if [ ! -s "$outfile" ]; then
  # file is empty or does not exist
    echo -e "$header" > "$outfile"
fi

for FILE in `cat its2_errfile_list.txt`; do
    # parse mapping rates from botwtie2 otuput
    jobid=$(echo $FILE | cut -d "_" -f1,2)
    sample=$(head -3 ../outfiles/$jobid'_its2_mapping_array_2024-07-24.out' | tail -1)
    con0=$(grep -oP '(?<=\s{4})\d+' $jobid'_its2_mapping_array_2024-07-24.err' | head -1)
    con1=$(grep -oP '(?<=\s{4})\d+' $jobid'_its2_mapping_array_2024-07-24.err' | head -2 | tail -1)
    con2=$(grep -oP '(?<=\s{4})\d+' $jobid'_its2_mapping_array_2024-07-24.err' | head -3 | tail -1)
    overall=$(grep -oP '\d+\.\d+%? overall alignment rate' $jobid'_its2_mapping_array_2024-07-24.err' | cut -d" " -f1)

    # Append data to output file
    echo -e "$sample\t$con0\t$con1\t$con2\t$overall" >> "$outfile"

done
```

 

``` bash
./parse_bowtie2_output_its2.sh
```

 

Conclusion: The read counts generated by the `countxpression` python
script are much higher than the mapping stats reported by Bowtie2. Upon
further investigation, it seems part of this is because the summary
Bowtie2 prints to the stdout actually reports counts of mapping read
PAIRS in the first several lines, rather than single reads (this is not
clear at all until you try to reconcile their numbers). However, the
read counts given by the python script are still too high despite this.
Looking back at the bowtie2 mapping command used (which was originally
modified from Veronica’s scripts), I noticed there is a “-k 5” argument
included that results in up to 5 “valid alignments” being reported per
read pair. When the -k argument is used, the supplemental alignments
(i.e. not the best alignment) are given a MAPQ score of 255. It also
appears bowtie reports unaligned reads in this mode with the YT:Z:UP
code. To remove these unnecessary lines from the alignment file, I ran
the following samtools command.

``` bash
for FILE in `ls *Itrimmed.bam`; do
     crun.samtools samtools view -@38 -h $FILE | grep -v "YT:Z:UP" | grep -vP "\s255\s" > ${FILE%.*}'_k1.sam'
     done
```

 

Then I reran the `countxpression_SB_advbioinf.py` and
`ParseExpression2BigTable_advbioinf.py` scripts on the new alignment
files. Note, because all of the alignments are going to be multiply
mapped, we have to first change the “columntoextract” variable in the
`ParseExpression2BigTable_advbioinf.py` script from “2” to “4”.
Otherwise, it will only summarize the counts from the 2nd column and all
entries will be 0.

``` bash
module load container_env
module load python2

cd $BASEDIR/its2_mapping

crun.python2 python2 ../scripts/countxpression_SB_advbioinf.py *k1.sam


crun.python2 python2 ../scripts/ParseExpression2BigTable_advbioinf.py ../../references/its2_strain_names.txt pver_pilot_its2_counts_merged_filtered.txt no_match *k1_counts.txt

    # Script usage:
      #sys.argv[1] Input file name #list of your "gene" names (list must have a header line!)
      #sys.argv[2] Output file name
      #sys.argv[3] Text to add when no match
      #sys.argv[4:] Any number of files to add columns from
```

 

Sum of columns of the `pver_pilot_its2_counts_merged_filtered.txt`
output file (mapped read totals for each sample):

    Sum of column 2: 371
    Sum of column 3: 408
    Sum of column 4: 150
    Sum of column 5: 48
    Sum of column 6: 420
    Sum of column 7: 84
    Sum of column 8: 118
    Sum of column 9: 144
    Sum of column 10: 278
    Sum of column 11: 106
    Sum of column 12: 110
    Sum of column 13: 272
    Sum of column 14: 82
    Sum of column 15: 245
    Sum of column 16: 32
    Sum of column 17: 63

 

Now I can recreate the ITS2 symbiont proportion plots:

``` r
library(tidyverse)
library(reshape2)
library(RColorBrewer)


setwd("C:/Users/jason/OneDrive/Documents/odu/projects/global_search_genomics/pver_gwas/")

# load count data from ParseExpression2BigTable_advbioinf.py
raw <- read_tsv("pver_pilot_its2_counts_merged_filtered.txt", col_names = TRUE)

# make column names shorter
col_names <- c("sym", "VATI_21", "VATI_22", "VATI_23", "VATI_24", "VATI_25", "VATI_26", "VATI_27", "VATI_28", "VATI_29", "VATI_30", "VATI_31", "VATI_32", "VATI_33", "VATI_34", "VATI_35", "VATI_36")

colnames(raw) <- col_names

# reshape data to long format for summarizing and plotting
counts <- melt(raw) %>% mutate(clade = substr(sym,1,1)) %>% rename(sample = variable, count = value)

# summarize counts by sample and clade
sum_table <- counts %>% group_by(sample, clade) %>% summarise(clade_count = sum(count))

# plot total counts
ggplot(sum_table, aes(fill=clade, y=clade_count, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Set1") +
  ylab("read count") +
  ggtitle("Raw Read Counts - ITS2 Symbiont Mapping - P. verrucosa - Vatia 2024") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_blank())

ggsave(plot = last_plot(), file = "2024_pver_pilot_symbiont_total_counts_plot.png", width = 8, height = 8, units = "in")

# calculate percent reads for each clade for each sample
total_reads <- counts %>% group_by(sample) %>% summarise(total_count = sum(count))

# merge data frames to add percent column to summary table
sum_table_perc <- left_join(sum_table, total_reads, by = "sample") %>% mutate(perc = clade_count/total_count)

# plot read counts as percentages
ggplot(sum_table_perc, aes(fill=clade, y=perc, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Set1") +
  ylab("% of mapped reads") +
  ggtitle("Read Proportions - ITS2 Symbiont Mapping - P. verrucosa - Vatia 2024") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_blank())


ggsave(plot = last_plot(), file = "2024_pver_pilot_symbiont_proportions_plot.png", width = 8, height = 8, units = "in")
```

 

![](C:/Users/jason/OneDrive/Documents/odu/projects/global_search_genomics/pver_gwas/2024_pver_pilot_symbiont_total_counts_plot.png)

![](C:/Users/jason/OneDrive/Documents/odu/projects/global_search_genomics/pver_gwas/2024_pver_pilot_symbiont_proportions_plot.png)
 

Convert .sam files to .bam files to save space:

``` bash
module load samtools

for FILE in *.sam; do
    BASENAME=`basename $FILE .sam`
    echo $BASENAME
    crun.samtools samtools view -Sb -@ 38 -O BAM -o $BASENAME.bam $FILE
done
```

 

## Create combined P. verrucosa + Cladocopium + Durusdinium reference

The references we will be using are the *Pocillopora verrucosa* genome
assembly published in Han et al. 2024 (GCA_020536085.1), the
*Cladocopium goreaui* genome assembly published in Chen et al. 2022
(accession GCA_947184155.1), and one of the *Durusdinium trenchii*
genome assemblies published in Dougan et al. 2022 (preprint, BioRxiv).  

Check to make sure P. verrucosa reference has 80 bp per line. It does.  

Add suffix to first entry of sequence name line for P. verrucosa
reference fasta to easily identify species:

``` bash
module load container_env
module load python2

crun python addsuffixtofastaseqnames.py Pverrucosa /cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF_036669915.1_ASM3666991v2_genomic.fna
```

 

Combine the references into one fasta file:

``` bash
cd /cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes

cat pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/GCF_036669915.1_ASM3666991v2_genom_suffixed.fasta GCA_947184155.1_Cgoreaui_SCF055-01_genom_suffixed.fasta Durusdinium_trenchii_SCF082.genome_rf.fa > combined_pver_cd_hologenome.fa
```

 

## Index reference file with Samtools (v1.19) and Bowtie2 (v2.4.1) and create sequence dictionary with Picard (v2.27.1)

``` bash
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
REFERENCE=$BASEDIR/references/genomes/combined_pver_cd_hologenome.fa   # This is the reference genome sequence to which we will map our reads
REFBASENAME="${REFERENCE%.*}" # Extract base file name from the reference for our output files

module load container_env
module load samtools

# samtools index
crun.samtools samtools faidx $REFERENCE

# create picard dictionary for later
module load picard/3.1.0
crun.picard picard CreateSequenceDictionary -R $REFERENCE -O $REFBASENAME'.dict'

# bowtie2 index
module load bowtie2
crun.bowtie2 bowtie2-build --threads 40 $REFERENCE $REFBASENAME
# takes about an 15 min to run depending on reference. One hour if only using one thread.
```

 

## Mapping to reference:

Run `hologenome_mapping_array.slurm`. This script will map each trimmed
fasta file to the hologenome reference file, convert each sam to bam,
index each bam, run idx stats, and then filter each bam for reads
mapping to Ahya contigs, creating a new bam file.

``` bash
#!/bin/bash

#SBATCH --job-name pver_hologenome_mapping_array_txgen_2024-08-02
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=28


## Load modules
module load container_env bowtie2

## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
FASTQDIR=$BASEDIR/pver_gwas_pilot/trimmed_fastq #path to trimmed fastq.gz files
OUTDIR=$BASEDIR/pver_gwas_pilot/bam
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_list_pver_pilot.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.
SAMPLETABLE=$BASEDIR/pver_gwas_pilot/sample_lists/fastq_table_pver_pilot.txt # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
FASTQ_SUFFIX_1=_f_paired_trim.fastq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.
FASTQ_SUFFIX_2=_r_paired_trim.fastq.gz # Suffix to raw fastq files. Use reverse reads with paired-end data.
REFBASENAME=combined_pver_cd_hologenome
#SCAFLIST=$BASEDIR/references/genomes/Ahyacinthus_scaffold_names_singleline.txt

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
SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2` # Define platform unit (PU), which is the lane number

echo $SAMPLE_UNIQ_ID

## Define the output path and file prefix
SAMPLEOUT=$OUTDIR/$SAMPLE_UNIQ_ID


## Run bowtie2 to map paired-end reads
# Note: we ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.
crun.bowtie2 bowtie2 -q --phred33 --very-sensitive -p 28 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA \
        -x $BASEDIR/references/genomes/$REFBASENAME -1 $FASTQDIR/$SAMPLE_UNIQ_ID$FASTQ_SUFFIX_1 -2 $FASTQDIR/$SAMPLE_UNIQ_ID$FASTQ_SUFFIX_2 \
        -S $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'

## Change modules
module unload bowtie2
module load container_env samtools

## Manipulate and summarise mapping with samtools
#sort and convert to BAM at the same time
crun.samtools samtools sort -@ 28 -o $SAMPLEOUT'_bt2_'$REFBASENAME'_sorted.bam' $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'

#remove SAM file
rm $SAMPLEOUT'_bt2_'$REFBASENAME'.sam'

#index sorted BAM
crun.samtools samtools index -@ 28 $SAMPLEOUT'_bt2_'$REFBASENAME'_sorted.bam'

#get summary stats
crun.samtools samtools idxstats -@ 28 $SAMPLEOUT'_bt2_'$REFBASENAME'_sorted.bam' > $SAMPLEOUT'_bt2_'$REFBASENAME'.idxstats'


# Extract mappings that are primary alignments only (no secondary or supplementary reads), with mapping score > 20, mapping length (CIGAR) > 20, and only on Ahya scaffolds
#crun.samtools samtools view -b -F 260 -q 20 -m 20 -@ 39 $SAMPLEOUT'_bt2_'$REFBASENAME'_sorted.bam' `cat $SCAFLIST`> $OUTDIR/ahya_bams/$SAMPLE_UNIQ_ID'_bt2_'$REFBASENAME'_sorted_primary_minq20_mlen20_ahya.bam'
```

 

Run the following script (`parse_bowtie2_output.sh`) to parse summary
mapping data from each file. The output file is called
`bowtie_mapping_summary.tsv`

``` bash
#!/bin/bash

ls *hologenome_mapping*.err > errfile_list.txt

# create data file and add header if it is empty
header="sample\tconcordantly_0_times\tconcordantly_1_time\tconcordantly_2_or_more_times\toverall_rate" #assign header value
outfile="bowtie_mapping_summary.tsv"

if [ ! -s "$outfile" ]; then
  # file is empty or does not exist
    echo -e "$header" > "$outfile"
fi

for FILE in `cat errfile_list.txt`; do
    # parse mapping rates from botwtie2 otuput
    jobid=$(echo $FILE | cut -d "_" -f1,2)
    sample=$(head -3 $jobid'_hologenome_mapping_array_pver_pilot_2024-08-02.out' | tail -1)
    con0=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pver_pilot_2024-08-02.err' | head -2 | tail -1)
    con1=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pver_pilot_2024-08-02.err' | head -3 | tail -1)
    con2=$(grep -oP '\d+\.\d+%?' $jobid'_hologenome_mapping_array_pver_pilot_2024-08-02.err' | head -4 | tail -1)
    overall=$(grep -oP '\d+\.\d+%? overall alignment rate' $jobid'_hologenome_mapping_array_pver_pilot_2024-08-02.err' | cut -d" " -f1)

    # Append data to output file
    echo -e "$sample\t$con0\t$con1\t$con2\t$overall" >> "$outfile"

done
```

 

``` bash
./parse_bowtie2_output.sh
```

 

    sample  concordantly_0_times    concordantly_1_time     concordantly_2_or_more_times    overall_rate
    2024_VATI_Pver_30_1_227H3WLT4-L008      26.56%  49.31%  24.13%  87.15%
    2024_VATI_Pver_31_1_227H3WLT4-L008      29.50%  46.39%  24.11%  86.23%
    2024_VATI_Pver_32_1_227H3WLT4-L008      30.38%  46.24%  23.38%  86.38%
    2024_VATI_Pver_33_1_227H3WLT4-L008      27.41%  49.35%  23.24%  87.21%
    2024_VATI_Pver_34_1_227H3WLT4-L008      29.85%  46.59%  23.56%  86.59%
    2024_VATI_Pver_35_1_227H3WLT4-L008      26.95%  49.37%  23.68%  87.30%
    2024_VATI_Pver_36_1_227H3WLT4-L008      25.53%  50.42%  24.04%  87.55%
    2024_VATI_Pver_21_1_227H3WLT4-L008      27.14%  47.87%  24.99%  86.66%
    2024_VATI_Pver_22_1_227H3WLT4-L008      27.76%  47.60%  24.64%  86.60%
    2024_VATI_Pver_23_1_227H3WLT4-L008      28.00%  47.36%  24.64%  86.45%
    2024_VATI_Pver_24_1_227H3WLT4-L008      28.95%  47.05%  24.00%  86.53%
    2024_VATI_Pver_25_1_227H3WLT4-L008      28.55%  47.84%  23.61%  86.49%
    2024_VATI_Pver_26_1_227H3WLT4-L008      27.77%  48.05%  24.18%  86.64%
    2024_VATI_Pver_27_1_227H3WLT4-L008      29.05%  47.11%  23.84%  86.53%
    2024_VATI_Pver_28_1_227H3WLT4-L008      31.60%  45.20%  23.19%  86.26%
    2024_VATI_Pver_29_1_227H3WLT4-L008      28.89%  47.25%  23.86%  85.95%

 

Move mapping summary table to `bam` directory:

``` bash
mv bowtie_mapping_summary.tsv ../bam/
```

 

## Merge files for the same individual (and same library!) that were sequenced in different batches

### Do NOT merge files for separately prepped libraries for the same individual.

We can skip this part since there are no replicates in this batch!  

|                                                                      |
|----------------------------------------------------------------------|
| \## CONTINUE WITH `GATK4_pipeline_jt.Rmd` FOR HOST SEQUENCE ANALYSIS |

 

## Separate out symbiont reads

### Create new bams for reads mapping to Cladocopium:

``` bash
mkdir $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/cgor_bams
```

 

`select_Cgor_alignments_array.slurm`

``` bash
#!/bin/bash

#SBATCH --job-name select_Cgor_alignments_array_2024-08-15
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
OUTDIR=$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/cgor_bams
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_bams_coordsorted_list.txt # Path to a bam list
REFBASENAME=PverCD
SCAFLIST=$BASEDIR/references/genomes/Cgoreaui_scaffold_names_singleline.txt



## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Define the output file name
CGOROUT=${SAMPLEFILE%combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam}$REFBASENAME'_dedup_primary_minq20_mlen20_cgor.bam'

# Load module
module load container_env samtools


#index sorted BAM
#crun.samtools samtools index -@ 20 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/$SAMPLEFILE

# Extract mappings that are primary alignments only (no unmapped or secondary/supplementary reads), with mapping score > 20, mapping length (CIGAR) > 20, and only on host scaffolds
crun.samtools samtools view -b -F 260 -q 20 -m 20 -@ 28 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/$SAMPLEFILE `cat $SCAFLIST`> $OUTDIR/$CGOROUT
```

 

Run samtools flagstat:

``` bash
for FILE in `ls *cgor.bam`; do
  echo ${FILE%227H*}
  crun.samtools samtools flagstat $FILE > ${FILE%.*}_flagstat.txt
done
```

 

Display total mapped and properly paired reads:

``` bash
for FILE in `ls *flagstat.txt`; do
  echo ${FILE%227H*}
  grep "0 mapped" $FILE
  grep "properly paired" $FILE
done
```

    2024_VATI_Pver_21_1_
    48133 + 0 mapped (100.00% : N/A)
    72 + 0 properly paired (0.15% : N/A)
    2024_VATI_Pver_22_1_
    54791 + 0 mapped (100.00% : N/A)
    128 + 0 properly paired (0.23% : N/A)
    2024_VATI_Pver_23_1_
    35734 + 0 mapped (100.00% : N/A)
    80 + 0 properly paired (0.22% : N/A)
    2024_VATI_Pver_24_1_
    56140 + 0 mapped (100.00% : N/A)
    54 + 0 properly paired (0.10% : N/A)
    2024_VATI_Pver_25_1_
    59798 + 0 mapped (100.00% : N/A)
    90 + 0 properly paired (0.15% : N/A)
    2024_VATI_Pver_26_1_
    62610 + 0 mapped (100.00% : N/A)
    48 + 0 properly paired (0.08% : N/A)
    2024_VATI_Pver_27_1_
    55138 + 0 mapped (100.00% : N/A)
    62 + 0 properly paired (0.11% : N/A)
    2024_VATI_Pver_28_1_
    44520 + 0 mapped (100.00% : N/A)
    52 + 0 properly paired (0.12% : N/A)
    2024_VATI_Pver_29_1_
    75717 + 0 mapped (100.00% : N/A)
    1804 + 0 properly paired (2.38% : N/A)
    2024_VATI_Pver_30_1_
    39624 + 0 mapped (100.00% : N/A)
    46 + 0 properly paired (0.12% : N/A)
    2024_VATI_Pver_31_1_
    43549 + 0 mapped (100.00% : N/A)
    60 + 0 properly paired (0.14% : N/A)
    2024_VATI_Pver_32_1_
    76736 + 0 mapped (100.00% : N/A)
    68 + 0 properly paired (0.09% : N/A)
    2024_VATI_Pver_33_1_
    45470 + 0 mapped (100.00% : N/A)
    42 + 0 properly paired (0.09% : N/A)
    2024_VATI_Pver_34_1_
    41965 + 0 mapped (100.00% : N/A)
    56 + 0 properly paired (0.13% : N/A)
    2024_VATI_Pver_35_1_
    50842 + 0 mapped (100.00% : N/A)
    28 + 0 properly paired (0.06% : N/A)
    2024_VATI_Pver_36_1_
    54614 + 0 mapped (100.00% : N/A)
    32 + 0 properly paired (0.06% : N/A)

 

Take a look at mapped sequences:

``` bash
crun.samtools samtools view 2024_VATI_Pver_21_1_227H3WLT4-L008_bt2_PverCD_dedup_primary_minq20_mlen20_cgor.bam | cut -f10 | less -S
```

    AAAGAAACAAGAAAGAAACAAGAAAGAAACAAGAAAAAAACAAGAAAGAAACAAGAAAAAAACAAGAAAGAAACAAGAAAGAAACAAGAAAGAAACAAAAAAGAAACAAGAAAGAAACAAGAAAGAAACAAGAAAGAAACAAGAAAGAAA
    GGTGCTGGTGCTGGTGCTGGTGCTGCTGGTGCTGCTGCTGGTGCTGGTGCTGCTGCTGCTGCTGGTGCTGCTGTTGCTGCTGGTGCTGCTGGTGCTGCTGGTGCTGCTGGTGCTGCTGCTGCTGCTGCTGGTGCTGCTGCTGCTGCTGCT
    GGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGAGAGGGGGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGAGGGGGTGGGGGGGTGGGGGGGGGGGAGGGGGGGAGGGGGGGGGGG
    AGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGGGGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGTGGGGGGGGGGGGGGGGTACGGGGGGGGGGGGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGAGAGGG
    CCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCC
    CCCCCCCCCCTCCCCCACCCCCCCCCCCCCCCCCCCCACCCCCCCTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    CCCTCTCTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    TCTCCCATGCTCGCAGGTCTTTGTCGGACTGCGCAATTGGCATACTTGCTGAAATTCCATACTTGCAGCTTCTCAATCTGGGCGTTCAGCGTGTTTTGGTCTGGCTCTGAAGTCACAGTCACGCTGTCTGATGCCTTCGCATGCATGGAG
    GCAAGTACGAAAAGGCTCAACACCATAGCTGGCATCACGAAACTTTCGGATATGCGTTGAACCTATGGAGCTGGCAGCAAAACACTGCGTCGGAAGTGACGCTTAAAAAATGCAAGGTCGGATATGAATGACATCATCCCAGATTTGCAC
    GGTATGTATGTATGTAGGTATGTATGTATGTGTGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATGTATG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTGTTGGGGGGGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGTTTC
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGGGGTGAAG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGAAGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGTGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAGGGGGGCGGTGGGGGGGGGGGGGGGGGGGGGGGGGGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGGGGCGGGGG
    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGGGGCGGGGGGGAGG

 

**Interpretation:** Upon inspection of the Cgor bam files, it seems most
of these alignments are improperly paired (e.g. forward read aligned to
a Cgor contig, but it’s mate is paired to a Pver contig). Looking at the
sequences, they seem to almost all be repeats as well, so I would not
trust any of these alignments to be true Cgor reads.  

### Create new bams for reads mapping to Durusdinium:

``` bash
mkdir $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/cgor_bams
```

 

`select_Dtre_alignments_array.slurm`

``` bash
#!/bin/bash

#SBATCH --job-name select_Dtre_alignments_array_2024-08-15
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
OUTDIR=$BASEDIR/pver_gwas_pilot/bam/dedup_bams2/dtre_bams
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_bams_coordsorted_list.txt # Path to a bam list
REFBASENAME=PverCD
SCAFBED=$BASEDIR/references/genomes/Dtrenchii.bed



## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Define the output file name
DTREOUT=${SAMPLEFILE%combined_pver_cd_hologenome_sorted_reheadered_qsorted_dedup_coordsorted.bam}$REFBASENAME'_dedup_primary_minq20_mlen20_dtre.bam'

# Load module
module load container_env samtools


#index sorted BAM
#crun.samtools samtools index -@ 20 $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/$SAMPLEFILE

# Extract mappings that are primary alignments only (no unmapped or secondary/supplementary reads), with mapping score > 20, mapping length (CIGAR) > 20, and only on host scaffolds
crun.samtools samtools view -b -F 260 -q 20 -m 20 -@ 28 --regions-file $SCAFBED $BASEDIR/pver_gwas_pilot/bam/dedup_bams2/$SAMPLEFILE > $OUTDIR/$DTREOUT
```

 

Upon inspection of the bam files, these alignments seem to be much
better than the Cgor alignments.

Run samtools flagstat:

``` bash
for FILE in `ls *dtre.bam`; do
  echo ${FILE%227H*}
  crun.samtools samtools flagstat $FILE > ${FILE%.*}_flagstat.txt
done
```

 

Display total mapped and properly paired reads:

``` bash
for FILE in `ls *flagstat.txt`; do
  echo ${FILE%227H*}
  grep "0 mapped" $FILE
  grep "properly paired" $FILE
done
```

    2024_VATI_Pver_21_1_
    420313 + 0 mapped (100.00% : N/A)
    367096 + 0 properly paired (87.34% : N/A)
    2024_VATI_Pver_22_1_
    485801 + 0 mapped (100.00% : N/A)
    427064 + 0 properly paired (87.91% : N/A)
    2024_VATI_Pver_23_1_
    246481 + 0 mapped (100.00% : N/A)
    216898 + 0 properly paired (88.00% : N/A)
    2024_VATI_Pver_24_1_
    159409 + 0 mapped (100.00% : N/A)
    138354 + 0 properly paired (86.79% : N/A)
    2024_VATI_Pver_25_1_
    528677 + 0 mapped (100.00% : N/A)
    457014 + 0 properly paired (86.44% : N/A)
    2024_VATI_Pver_26_1_
    239005 + 0 mapped (100.00% : N/A)
    208912 + 0 properly paired (87.41% : N/A)
    2024_VATI_Pver_27_1_
    164978 + 0 mapped (100.00% : N/A)
    140202 + 0 properly paired (84.98% : N/A)
    2024_VATI_Pver_28_1_
    153756 + 0 mapped (100.00% : N/A)
    127434 + 0 properly paired (82.88% : N/A)
    2024_VATI_Pver_29_1_
    276251 + 0 mapped (100.00% : N/A)
    235130 + 0 properly paired (85.11% : N/A)
    2024_VATI_Pver_30_1_
    127109 + 0 mapped (100.00% : N/A)
    108742 + 0 properly paired (85.55% : N/A)
    2024_VATI_Pver_31_1_
    224169 + 0 mapped (100.00% : N/A)
    191534 + 0 properly paired (85.44% : N/A)
    2024_VATI_Pver_32_1_
    319386 + 0 mapped (100.00% : N/A)
    269216 + 0 properly paired (84.29% : N/A)
    2024_VATI_Pver_33_1_
    177849 + 0 mapped (100.00% : N/A)
    147418 + 0 properly paired (82.89% : N/A)
    2024_VATI_Pver_34_1_
    256496 + 0 mapped (100.00% : N/A)
    216330 + 0 properly paired (84.34% : N/A)
    2024_VATI_Pver_35_1_
    92904 + 0 mapped (100.00% : N/A)
    78092 + 0 properly paired (84.06% : N/A)
    2024_VATI_Pver_36_1_
    121881 + 0 mapped (100.00% : N/A)
    105248 + 0 properly paired (86.35% : N/A)

## Deduplicate (Picard v3.1.0) and clip overlapping read pairs (BamUtil v1.0.15)

Make new bam list to be used for deduplication:

``` bash
cd $BASEDIR/bam

ls *sorted.bam > ../sample_lists/bam_list_for_dedup.txt
```

 

Script for all (merged and single) bams (dedup_clip_array.slurm):

``` bash
#!/bin/bash
#SBATCH --job-name dedup_clip_array_2024-08-06
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=28


## Load modules
module load container_env picard
# BamUtil is installed locally and doesn't need a module/container to run. Execute with "bamutil" command

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/bam_list_for_dedup.txt
REFNAME=pvercd-hologenome # Reference name to add to output files
PICARD='crun.picard picard'
BAMUTIL=~/.local/bin/bamutil #bamutil -h for usage (note: this executable was changed to "bamutil" from the original "bam" for clarity)

# Make directory for output
mkdir $BASEDIR/pver_gwas_pilot/bam/dedup_bams
mkdir $BASEDIR/pver_gwas_pilot/bam/dedup_bams/clipped_bams

## Loop over each sample
# for SAMPLEBAM in `cat $BAMLIST`; do

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
echo Sample bam is $SAMPLEBAM


## Remove duplicates and print dupstat file
$PICARD -Xmx300g MarkDuplicates \
  I=$BASEDIR'/pver_gwas_pilot/bam/'$SAMPLEBAM \
  O=$BASEDIR'/pver_gwas_pilot/bam/dedup_bams/'${SAMPLEBAM%.*}'_dedup.bam' \
  M=$BASEDIR'/pver_gwas_pilot/bam/dedup_bams/'${SAMPLEBAM%.*}'_dupstat.txt' \
  VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

## Clip overlapping paired end reads (only necessary for paired-end data)
$BAMUTIL clipOverlap --in $BASEDIR'/pver_gwas_pilot/bam/dedup_bams/'${SAMPLEBAM%.*}'_dedup.bam' --out $BASEDIR'/pver_gwas_pilot/bam/dedup_bams/clipped_bams/'${SAMPLEBAM%.*}'_dedup_overlapclipped.bam' --stats

echo 'done-zo woot!'
```

 

``` bash
sbatch dedup_clip_array.slurm
```

 

Check summary of duplicates for each file:

``` bash
for FILE in `ls *dupstat.txt`; do grep -A 1 "LIBRARY" $FILE | grep -v "LIBRARY"; done
```

 

    LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED     SECONDARY_OR_SUPPLEMENTARY_RDS  UNMAPPED_READS  UNPAIRED_READ_DUPLICATES        READ_PAIR_DUPLICATES   READ_PAIR_OPTICAL_DUPLICATES    PERCENT_DUPLICATION     ESTIMATED_LIBRARY_SIZE
    2024_VATI_Pver_21       36983145        247006527       0       81732269        27300978        83411216        48188053    0.365583        492679476
    2024_VATI_Pver_22       39157725        251995561       0       84018949        29179454        92905508        59852342    0.395822        492404623
    2024_VATI_Pver_23       33533776        221112947       0       74542518        24442737        79222453        44257118    0.384412        386140390
    2024_VATI_Pver_24       34871082        219228232       0       73683510        25549836        87918884        53410027    0.425472        340948168
    2024_VATI_Pver_25       27898232        182869676       0       61503878        19642829        72101489        47099547    0.416235        321840390
    2024_VATI_Pver_26       25243806        164779816       0       54699392        17219839        60689897        42566124    0.390638        370227369
    2024_VATI_Pver_27       27002529        171467373       0       57583329        18786074        64092784        40914971    0.397288        322734019
    2024_VATI_Pver_28       29843587        179911952       0       62050905        21093615        66759701        42597277    0.396782        342921670
    2024_VATI_Pver_29       34713874        223747812       0       78825490        25843195        92002388        51535021    0.43518         306446353
    2024_VATI_Pver_30       23769838        165072990       0       52166900        16491913        70826100        45280737    0.446841        239337862
    2024_VATI_Pver_31       30640894        201052214       0       69114360        22393057        78917304        43480119    0.416475        295552164
    2024_VATI_Pver_32       38795915        239386675       0       81600393        29146185        94692537        54647775    0.422226        362023035
    2024_VATI_Pver_33       16929501        109617561       0       34633729        10541017        40768208        33065480    0.389887        354429149
    2024_VATI_Pver_34       34994577        220201738       0       73606453        25958842        90567770        50197785    0.435623        298706381
    2024_VATI_Pver_35       15433097        103782230       0       32428883        9528990         37813445        28625399    0.381869        281786784
    2024_VATI_Pver_36       16296069        112517752       0       34315793        10387839        44627299        33250655    0.412886        249039936

 

Create new list for deduped and overlap clipped bams:

``` bash
cd $BASEDIR/pver_gwas_pilot/bam/dedup_bams/clipped_bams/

ls *.bam > $BASEDIR/pver_gwas_pilot/sample_lists/dedup_clipped_bam.list
ls /cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/bam/dedup_bams/clipped_bams/*overlapclipped.bam > /cm/shared/courses/dbarshis/barshislab/jtoy/pver_gwas_pilot/sample_lists/dedup_clipped_bam_fullpath.list # For use in the indel realignment step, this list must have the ".list" extension and include full paths
```

 

## Indel Realignment using GATK (v3.8)

Unlike other variant detector programs (like the GATK Haplotype Caller
or Freebayes), ANGSD does not realign reads during its analysis. Because
it can be difficult to distinguish indels from SNPs at the end of reads
if each alignment is considered separately, indels may interfere with
genotype likelihood estimation. We therefore recommend running your bam
files through a program that realigns reads around indels prior to
running ANGSD. The GATK IndelRealigner takes all the aligned sequences
from all samples into account to validate the indels discovered from the
mapping process and then realigns each read locally. This is not a
mandatory step and tends to be very time-consuming if you have a large
dataset, but the code is provided here.  

First index the dedupped/clipped BAMs:
`index_dedup-clipped_bams_array.slurm`

``` bash
#!/bin/bash

#SBATCH --job-name index_dedup-clipped_bams_array_2024-08-07
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=38



BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_clipped_bam.list # Path to a list of merged, deduplicated, and overlap clipped bam files.


## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE


## Index BAM
# Load samtools module
module load container_env samtools

#index sorted BAM
crun.samtools samtools index -@ 38 $BASEDIR/pver_gwas_pilot/bam/dedup_bams/clipped_bams/$SAMPLEFILE

echo 'done-zo!'
```

``` bash
sbatch index_dedup-clipped_bams_array.slurm
```

 

Realign reads around indels using tools from GATK (v3.8-1-0-gf15c1c3ef)
This step gets run across all samples at once using a single command so
can’t be split into an array job. `indel_realignment.slurm`

``` bash
#!/bin/bash

#SBATCH --job-name indel_realignment_2024-08-07
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=38

BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/dedup_clipped_bam_fullpath.list # Path to a list of merged, deduplicated, and overlap clipped bam files. Must be indexed. Full paths should be included. This file has to have a suffix of ".list".
REFERENCE=$BASEDIR/references/genomes/combined_pver_cd_hologenome.fa


## Realign around in-dels
# This is done across all samples at once!
module load container_env gatk/3.8

## Create list of potential in-dels
crun.gatk gatk3 -Xmx200g \
-T RealignerTargetCreator \
-R $REFERENCE \
-I $BAMLIST \
-o $BASEDIR'/pver_gwas_pilot/bam/dedup_bams/clipped_bams/all_samples_for_indel_realigner.intervals' \
-drf BadMate \
-nt 38 # number of threads to use

## Run the indel realigner tool
crun.gatk gatk3 -Xmx100g \
-T IndelRealigner \
-R $REFERENCE \
-I $BAMLIST \
-targetIntervals $BASEDIR'/pver_gwas_pilot/bam/dedup_bams/clipped_bams/all_samples_for_indel_realigner.intervals' \
--consensusDeterminationModel USE_READS  \
--nWayOut _realigned.bam
```

 

Move realigned bams to new directory:

``` bash
mkdir realigned_bams

mv *realigned.bam realigned_bams
```

Make list of realigned bams:

``` bash
cd $BASEDIR/pver_gwas_pilot/bam/dedup_bams/clipped_bams/realigned_bams

ls *realigned.bam > $BASEDIR/pver_gwas_pilot/sample_lists/realigned_bam_list.txt
```

# Filter BAMs to only include Pver alignments

Make scaffold list for host:

``` bash
cd /cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes

grep -P '^>.+Pverrucosa' combined_pver_cd_hologenome.fa | sed 's/^>//' > Pver_scaffold_names.txt
grep -P '^>.+Pverrucosa' combined_pver_cd_hologenome.fa | sed 's/^>//' | tr '\n' ' ' > Pver_scaffold_names_singleline.txt
```

 

Filter bams using `select_host_alignments_array.slurm` script:

``` bash
#!/bin/bash

#SBATCH --job-name select_host_alignments_array_2024-08-07
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-16%16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time 5-00:00:00
#SBATCH --cpus-per-task=28



## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy
OUTDIR=$BASEDIR/pver_gwas_pilot/bam/dedup_bams/clipped_bams/pver_bams
SAMPLELIST=$BASEDIR/pver_gwas_pilot/sample_lists/realigned_bam_list.txt # Path to a bam list
REFBASENAME=PverCD
SCAFLIST=$BASEDIR/references/genomes/Pver_scaffold_names_singleline.txt


## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Define the output path and file prefix
HOSTOUT=${SAMPLEFILE%combined_pver_cd_hologenome_sorted_dedup_overlapclipped.bam}$REFBASENAME'_primary_minq20_mlen20_pver.bam'

# Load module
module load container_env samtools


#index sorted BAM
crun.samtools samtools index -@ 28 $BASEDIR/pver_gwas_pilot/bam/dedup_bams/clipped_bams/$SAMPLEFILE

# Extract mappings that are primary alignments only (no unmapped or secondary/supplementary reads), with mapping score > 20, mapping length (CIGAR) > 20, and only on host scaffolds
crun.samtools samtools view -b -F 260 -q 20 -m 20 -@ 28 $BASEDIR/pver_gwas_pilot/bam/dedup_bams/clipped_bams/$SAMPLEFILE `cat $SCAFLIST`> $OUTDIR/$HOSTOUT
```

``` bash
sbatch select_host_alignments_array.slurm
```

 

Create a new list of host-only BAM files:

``` bash
ls *mlen20_pver.bam > $BASEDIR/pver_gwas_pilot/sample_lists/pver_bam_list.txt
```

 

Count how many reads are left for each sample:

``` bash
for file in `cat pver_bam_list`; do
  crun.samtools samtools view $file | wc -l
done
```

 

## Estimate read depth of .bam files using samtools (v1.14)

 

First, run samtools depth to get depth per sample per position.

``` bash
BAMLIST=$BASEDIR/pver_gwas_pilot/sample_lists/pver_bam_list.txt

cd $BASEDIR/pver_gwas_pilot/bam/dedup_bams/clipped_bams/pver_bams

for SAMPLEBAM in `cat $BAMLIST`; do
    echo $SAMPLEBAM
    ## Count per position depth per sample (-aa = Output absolutely all positions, including zero depth positions and unused ref seqs)
    crun.samtools samtools depth -aa -@ 40 $SAMPLEBAM | gzip > ${SAMPLEBAM%.*}'.depth.gz'
done

mkdir samtools_depth_out
mv $BASEDIR/bam/realigned_bams/*.depth.gz samtools_depth_out
```

 

Download and check data from TxGen
================
Jason Toy
2024-05-28

## Extract original links from email

Open raw html code for email sent from TxGen. Use Notepad++ to isolate
original download links (not SafeLinks) for each file. Then add wget
commands to the start of each line. Run these commands to download the
files:

``` bash
wget https://download.txgen.tamu.edu/dan.barshis/240430_LH00550_0004_B22JWTTLT3_23313Brs_Join/240430_LH00550_0004_B22JWTTLT3-23313Brs_fastq.tar
wget https://download.txgen.tamu.edu/dan.barshis/240430_LH00550_0004_B22JWTTLT3_23313Brs_Join/240430_LH00550_0004_B22JWTTLT3-23313Brs_fastq.tar.md5
wget https://download.txgen.tamu.edu/dan.barshis/240430_LH00550_0004_B22JWTTLT3_23313Brs_Join/240430_LH00550_0004_B22JWTTLT3-23313Brs_fastqc.tar
wget https://download.txgen.tamu.edu/dan.barshis/240430_LH00550_0004_B22JWTTLT3_23313Brs_Join/240430_LH00550_0004_B22JWTTLT3-23313Brs_fastqc.tar.md5
wget https://download.txgen.tamu.edu/dan.barshis/240430_LH00550_0004_B22JWTTLT3_23313Brs_Join/240430_LH00550_0004_B22JWTTLT3-23313Brs_samplesheets.tar
wget https://download.txgen.tamu.edu/dan.barshis/240430_LH00550_0004_B22JWTTLT3_23313Brs_Join/240430_LH00550_0004_B22JWTTLT3-23313Brs_samplesheets.tar.md5
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0005_A22JY72LT3_23313Brs_Join/240507_LH00550_0005_A22JY72LT3-23313Brs_fastq.tar
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0005_A22JY72LT3_23313Brs_Join/240507_LH00550_0005_A22JY72LT3-23313Brs_fastq.tar.md5
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0005_A22JY72LT3_23313Brs_Join/240507_LH00550_0005_A22JY72LT3-23313Brs_fastqc.tar
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0005_A22JY72LT3_23313Brs_Join/240507_LH00550_0005_A22JY72LT3-23313Brs_fastqc.tar.md5
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0005_A22JY72LT3_23313Brs_Join/240507_LH00550_0005_A22JY72LT3-23313Brs_samplesheets.tar
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0005_A22JY72LT3_23313Brs_Join/240507_LH00550_0005_A22JY72LT3-23313Brs_samplesheets.tar.md5
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0006_B22JY53LT3_23313Brs_Join/240507_LH00550_0006_B22JY53LT3-23313Brs_fastq.tar
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0006_B22JY53LT3_23313Brs_Join/240507_LH00550_0006_B22JY53LT3-23313Brs_fastq.tar.md5
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0006_B22JY53LT3_23313Brs_Join/240507_LH00550_0006_B22JY53LT3-23313Brs_fastqc.tar
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0006_B22JY53LT3_23313Brs_Join/240507_LH00550_0006_B22JY53LT3-23313Brs_fastqc.tar.md5
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0006_B22JY53LT3_23313Brs_Join/240507_LH00550_0006_B22JY53LT3-23313Brs_samplesheets.tar
wget https://download.txgen.tamu.edu/dan.barshis/240507_LH00550_0006_B22JY53LT3_23313Brs_Join/240507_LH00550_0006_B22JY53LT3-23313Brs_samplesheets.tar.md5
```

 

## Make md5sum.txt file

MD5 sums for each file are provided, but the `md5sum.txt` file that is
required to run the md5sum check on all files is not, so we need to
create a new file in the format of `md5sum_code    original_file_name`:

``` bash
for FILE in `ls *.md5`; do
    echo -e `cat $FILE`'\t'${FILE%.*} >> md5sum.txt
done
```

 

## Run md5sum check

Now we can run the check on the files with the `md5sum` command:

``` bash
md5sum -c md5sum.txt
```

    240430_LH00550_0004_B22JWTTLT3-23313Brs_fastqc.tar: OK
    240430_LH00550_0004_B22JWTTLT3-23313Brs_fastq.tar: OK
    240430_LH00550_0004_B22JWTTLT3-23313Brs_samplesheets.tar: OK
    240507_LH00550_0005_A22JY72LT3-23313Brs_fastqc.tar: OK
    240507_LH00550_0005_A22JY72LT3-23313Brs_fastq.tar: OK
    240507_LH00550_0005_A22JY72LT3-23313Brs_samplesheets.tar: OK
    240507_LH00550_0006_B22JY53LT3-23313Brs_fastqc.tar: OK
    240507_LH00550_0006_B22JY53LT3-23313Brs_fastq.tar: OK
    240507_LH00550_0006_B22JY53LT3-23313Brs_samplesheets.tar: OK

 

## Unpack all tar archives

To get gzipped fastq files

    for FILE in `ls *.tar`; do
      tar -xvf $FILE
    done

 

## Run FastQC on each fastq file by running array script `fastqc_array.slurm`

``` bash
#!/bin/bash

#SBATCH --job-name fastqc_array_2024-06-03
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-830%24
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time 5-00:00:00


## Load modules
module load container_env
module load fastqc/0.11.9

## Define some variables
BASEDIR=/cm/shared/courses/dbarshis/barshislab/jtoy/
RAWDATA=$BASEDIR/raw_sequence_data/ahya_gwas_2024-05-22_txgen #path to raw fq.gz files
OUTDIR=$RAWDATA/fastqc_jt/
SAMPLELIST=$RAWDATA/fastq_list.txt #Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.

## Keep a record of the Job ID
echo $SLURM_JOB_ID

## Select the SAMPLE from the SAMPLELIST
SAMPLEFILE=`head $SAMPLELIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## Keep record of sample file
echo $SAMPLEFILE

## Run FastQC
crun.fastqc fastqc -t 30 -o $OUTDIR $SAMPLEFILE
```

 

## Compile FastQC outputs into one report with MultiQC

``` bash
salloc --partition=main --nodes=1 --exclusive

module load container_env multiqc

cd /cm/shared/courses/dbarshis/barshislab/jtoy/raw_sequence_data/ahya_gwas_2024-05-22_txgen/fastqc_jt

crun.multiqc multiqc --interactive --filename multiqc_report_ahya_txgen_2024-05-22_fastqc_interactive_jt .
# --interactive forces the creation of an interactive html file even when sample sizes are high (instead of a flat file)
```

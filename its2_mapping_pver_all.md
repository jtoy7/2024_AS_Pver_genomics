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



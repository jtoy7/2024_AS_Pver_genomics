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

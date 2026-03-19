Calculation of standard population genetics metrics using VCFtools
================
**AUTHOR:** Jason A. Toy  
**DATE:** 2026-03-19 <br><br>


## Prepare dataset

<br>

For these analyses, we will be starting with the `pver_all_QDPSBfiltered_genotypes.vcf.gz` VCF file in the `/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf` directory, which has all original samples and is only filtered for quality, depth, and strand bias (bcftools filters: -e 'QUAL < 30 || INFO/MQ < 40 || INFO/DP < 792 || INFO/DP > 25286 || INFO/QD < 2.0 || INFO/FS > 60.0 || INFO/SOR > 3.0').

<br>

Start by creating a single-column sample list of clone-pruned P. acuta samples from the existing keep file. Also remove any "Extra" samples:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf

cut -f2 keep_samples_Pacuta_only.txt > keep_samples_Pacuta_only.samples

grep -v 'X' keep_samples_Pacuta_only.samples > keep_samples_Pacutaonly_noExtras.samples
```
This leaves 135 samples.

<br>

Then subset from the QDPSBfiltered VCF and keep only biallelic SNPs:
```bash
module load bcftools

cd /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf

crun.bcftools bcftools view \
  -S keep_samples_Pacutaonly_noExtras.samples \
  -v snps -m2 -M2 \
  -Ou \
  pver_all_QDPSBfiltered_genotypes.vcf.gz | \
bcftools +fill-tags -Ou -- -t AC,AN | \
bcftools view \
  -e 'AC=0 || AC=AN' \
  -Oz -o pver_Pacuta_clonepruned_snps_preMISS.vcf.gz

bcftools index pver_Pacuta_clonepruned_snps_preMISS.vcf.gz
```

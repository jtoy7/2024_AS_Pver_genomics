Calculation of standard population genetics metrics using pixy
================
**AUTHOR:** Jason A. Toy  
**DATE:** 2026-03-23 <br><br>


## Prepare dataset

<br>

`pixy` requires a VCF that includes invariant sites in addition to variant sites (this is how it accounts for missing data), so unfortunatly, we have to go back and recall genotypes from gVCFs in `GATK` using `GenotypeGVCFs` and the `--include-non-variant-sites` or `-all-sites` flag.

<br>

### Create sample/population files
Start by creating a few sample/population files we'll need. First up, create a single-column sample list of clone-pruned P. acuta samples from the existing keep file. Also remove any "Extra" samples. This file is for filtering with vcftools/bcftools and will be called 'keep_samples_Pacutaonly_noExtras.samples':
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf

cut -f2 keep_samples_Pacuta_only.txt > keep_samples_Pacuta_only.samples

grep -v 'X' keep_samples_Pacuta_only.samples > keep_samples_Pacutaonly_noExtras.samples
```
This leaves 135 samples.

<br>

Next, we'll create a few files to use with pixy:
`pixy_all.pop.tsv` for whole-dataset π
`pixy_location.pop.tsv` for per-location π and all pairwise location FST/dxy
`pixy_island.pop.tsv` for Ofu vs. Tutuila π/FST/dxy
```bash
# whole dataset = one population called ALL
awk 'BEGIN{OFS="\t"} {print $2, "ALL"}' keep_samples_Pacuta_only.txt | grep -v 'X' > pixy_all.pop.tsv

# per-location populations
awk 'BEGIN{OFS="\t"} {split($2,a,"_"); print $2, a[2]}' keep_samples_Pacuta_only.txt | grep -v 'X' > pixy_location.pop.tsv

# island populations
awk 'BEGIN{OFS="\t"} {
  split($2,a,"_");
  island = (a[2]=="OFU3" || a[2]=="OFU6") ? "Ofu" : "Tutuila";
  print $2, island
}' keep_samples_Pacuta_only.txt | grep -v 'X' > pixy_island.pop.tsv

# ternary operator is essentially a compact if/else statement: "condition ? value_if_true : value_if_false"
# so if Location is OFU3 or OFU6 assign "Ofu", otherwise assign "Tutuila"
```

### Regenerate all-sites VCFs from existing GenomicsDBs and gVCFs
Create and run a modified `GenotypeGVCFs_array.slurm` script using the `--include-non-variant-sites` or `-all-sites` flag and the `-L` region flag (apparently necessary to prevent issues when calling all sites).
```bash
#!/bin/bash

#SBATCH --job-name=GenotypeGVCFs_pver_allsites_2026-03-23
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
OUTDIR=$BASEDIR/pver_gwas/hologenome_mapped_all/vcf/allsites_vcf
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
  --include-non-variant-sites \
  -L ${SCAF}
  -O $OUTDIR/${SCAF}_genotypes.vcf.gz
```
The longest scaffold took about 24 hrs to complete.

<br>

Check output to make sure invariant sites are included:
```bash
module load bcftools

crun.bcftools bcftools view -H /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf_allsites/NW_027078162.1_Pverrucosa_allsites_genotypes.vcf.gz | head | less -S
```
```
NW_027078162.1_Pverrucosa       1       .       G       .       .       .       .       GT      ./.     ./.     ./.     ./.     ./. >
NW_027078162.1_Pverrucosa       2       .       G       .       .       .       .       GT      ./.     ./.     ./.     ./.     ./. >
NW_027078162.1_Pverrucosa       3       .       G       .       .       .       .       GT      ./.     ./.     ./.     ./.     ./. >
NW_027078162.1_Pverrucosa       4       .       C       .       .       .       .       GT      ./.     ./.     ./.     ./.     ./. >
NW_027078162.1_Pverrucosa       5       .       G       .       .       .       .       GT      ./.     ./.     ./.     ./.     ./. >
```
Shows that non-variant positions were emitted.

<br>

```bash
crun.bcftools bcftools view -H /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf_allsites/NW_027078162.1_Pverrucosa_allsites_genotypes.vcf.gz | grep '0/0' | head | less -S
```
```
.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. 0/0:2:6 ./.:.:. ./.:.:. >
.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. 0/0:2:6 ./.:.:. ./.:.:. >
.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. 0/0:2:6 ./.:.:. ./.:.:. >
.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. 0/0:2:6 ./.:.:. ./.:.:. >
.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. ./.:.:. 0/0:2:6 ./.:.:. ./.:.:. >
```
Shows that non-variant positions were actually called.

<br>

### Filter the all-sites VCF for each scaffold

Variant and non-variant sites will need to be filtered separately and then recombined.

- Start with hard filters for quality, depth, and strand bias.
- Then subset to the clone-pruned, P. acuta-only sample set
- Then apply missingness filter (max missingness per site = 20%)

Run this as an array job applying filters to each scaffold separately:
```bash
#!/bin/bash

#SBATCH --job-name=filter_allsites_vcf_for_pixy_2026-03-24
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-52%30
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=10

set -euo pipefail

module load bcftools vcftools

BASEDIR=/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all
INDIR=$BASEDIR/vcf_allsites
OUTDIR=$BASEDIR/vcf_allsites/filtered
SAMPLES=$BASEDIR/vcf/keep_samples_Pacutaonly_noExtras.samples
SCAFLIST=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/genome_regions.list

mkdir -p "$OUTDIR"

TASK_ID=${SLURM_ARRAY_TASK_ID}
SCAF=$(sed -n "${TASK_ID}p" "$SCAFLIST")
INVCF="${INDIR}/${SCAF}_allsites_genotypes.vcf.gz"

echo "Processing scaffold: $SCAF"
echo "Input VCF: $INVCF"
echo "Sample list: $SAMPLES"

if [[ ! -f "$INVCF" ]]; then
  echo "ERROR: Input VCF not found: $INVCF" >&2
  exit 1
fi

# 1) Filter for variant SNPs using the same hard filters from variant-only workflow,
#    then subset to clone-pruned, Pacuta-only samples.
crun.bcftools bcftools view \
  -v snps -m2 -M2 \
  "$INVCF" -Ou | \
crun.bcftools bcftools filter \
  -e 'QUAL < 30 || INFO/MQ < 40 || INFO/DP < 792 || INFO/DP > 25286 || INFO/QD < 2.0 || INFO/FS > 60.0 || INFO/SOR > 3.0' \
  -Ou | \
crun.bcftools bcftools view --threads 10 \
  -S "$SAMPLES" \
  -Oz -o "${OUTDIR}/${SCAF}.variant.subset.vcf.gz"

crun.bcftools bcftools index "${OUTDIR}/${SCAF}.variant.subset.vcf.gz"

# 2) Pull invariant sites, then subset to clone-pruned Pacuta-only samples.
#    No QUAL/QD/FS/SOR filtering here.
crun.bcftools bcftools view \
  -v none \
  "$INVCF" -Ou | \
crun.bcftools bcftools filter \
  -e 'INFO/DP < 792 || INFO/DP > 25286' \
  -Ou | \
crun.bcftools bcftools view --threads 10 \
  -S "$SAMPLES" \
  -Oz -o "${OUTDIR}/${SCAF}.invariant.subset.vcf.gz"

crun.bcftools bcftools index "${OUTDIR}/${SCAF}.invariant.subset.vcf.gz"

# 3) Apply the same missingness filter to each subset separately.
crun.vcftools vcftools \
  --gzvcf "${OUTDIR}/${SCAF}.variant.subset.vcf.gz" \
  --max-missing 0.8 \
  --recode --stdout | \
crun.bcftools bgzip -c > "${OUTDIR}/${SCAF}.variant.subset.miss80.vcf.gz"

crun.vcftools vcftools \
  --gzvcf "${OUTDIR}/${SCAF}.invariant.subset.vcf.gz" \
  --max-missing 0.8 \
  --recode --stdout | \
crun.bcftools bgzip -c > "${OUTDIR}/${SCAF}.invariant.subset.miss80.vcf.gz"

crun.bcftools bcftools index "${OUTDIR}/${SCAF}.variant.subset.miss80.vcf.gz"
crun.bcftools bcftools index "${OUTDIR}/${SCAF}.invariant.subset.miss80.vcf.gz"

# 4) Recombine the filtered variant and invariant records into a final pixy-ready all-sites VCF.
crun.bcftools bcftools concat --threads 10 \
  --allow-overlaps \
  "${OUTDIR}/${SCAF}.variant.subset.miss80.vcf.gz" \
  "${OUTDIR}/${SCAF}.invariant.subset.miss80.vcf.gz" \
  -Ou | \
crun.bcftools bcftools sort --threads 10 \
  -Oz -o "${OUTDIR}/${SCAF}.pixy_ready.vcf.gz"

crun.bcftools bcftools index "${OUTDIR}/${SCAF}.pixy_ready.vcf.gz"

echo "Finished scaffold: $SCAF"

```
Runtime:

<br>

Check how many nonvariant sites were removed by filters:
```bash
# variant sites
bcftools view -H -v snps -m2 -M2 input.vcf.gz | wc -l; \
bcftools view -v snps -m2 -M2 input.vcf.gz -Ou | bcftools filter -e 'INFO/DP < 792 || INFO/DP > 25286' -Ou | bcftools view -H | wc -l

# invariant sites
bcftools view -H -v none input.vcf.gz | wc -l; \
bcftools view -v none input.vcf.gz -Ou | bcftools filter -e 'INFO/DP < 792 || INFO/DP > 25286' -Ou | bcftools view -H | wc -l
```
As for-loop across all scaffolds:
```bash
SCAFLIST=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/genome_regions.list
INDIR=/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf_allsites

echo -e "scaffold\tvar_before\tvar_after\tvar_pct_removed\tinv_before\tinv_after\tinv_pct_removed"

while read SCAF; do
  INVCF="${INDIR}/${SCAF}_allsites_genotypes.vcf.gz"

  vb=$(bcftools view -H -v snps -m2 -M2 "$INVCF" | wc -l)
  va=$(bcftools view -v snps -m2 -M2 "$INVCF" -Ou | \
       bcftools filter -e 'INFO/DP < 792 || INFO/DP > 25286' -Ou | \
       bcftools view -H | wc -l)

  ib=$(bcftools view -H -v none "$INVCF" | wc -l)
  ia=$(bcftools view -v none "$INVCF" -Ou | \
       bcftools filter -e 'INFO/DP < 792 || INFO/DP > 25286' -Ou | \
       bcftools view -H | wc -l)

  vp=$(awk -v b=$vb -v a=$va 'BEGIN{if(b>0) printf "%.2f", (b-a)/b*100; else print "NA"}')
  ip=$(awk -v b=$ib -v a=$ia 'BEGIN{if(b>0) printf "%.2f", (b-a)/b*100; else print "NA"}')

  echo -e "${SCAF}\t${vb}\t${va}\t${vp}\t${ib}\t${ia}\t${ip}"

done < "$SCAFLIST"
```

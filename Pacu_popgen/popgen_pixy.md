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
The longest scaffold took 40.5 hrs to complete.

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

Run this as an array job applying filters to each scaffold separately. Note that VCFtools is installed in a module called `dosage_convertor`.
<br>

`filter_allsites_vcf_for_pixy_array.slurm`:
```bash
#!/bin/bash

#SBATCH --job-name=filter_allsites_vcf_for_pixy_array_2026-03-25
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

module load bcftools

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
  -S "$SAMPLES" -a \
  -Oz -o "${OUTDIR}/${SCAF}.variant.subset.vcf.gz"

crun.bcftools bcftools index "${OUTDIR}/${SCAF}.variant.subset.vcf.gz"

# 2) Pull invariant sites, then subset to clone-pruned Pacuta-only samples.
#    No QUAL/QD/FS/SOR filtering here; just DP filtering.
crun.bcftools bcftools view \
  -i 'ALT="."' \
  "$INVCF" -Ou | \
crun.bcftools bcftools filter \
  -e 'INFO/DP < 792 || INFO/DP > 25286' \
  -Ou | \
crun.bcftools bcftools view --threads 10 \
  -S "$SAMPLES" \
  -Oz -o "${OUTDIR}/${SCAF}.invariant.subset.vcf.gz"

crun.bcftools bcftools index "${OUTDIR}/${SCAF}.invariant.subset.vcf.gz"

# 3) Apply the same missingness filter to each subset separately.
crun.bcftools bcftools view --threads 10 \
  -i 'F_MISSING<=0.2' \
  "${OUTDIR}/${SCAF}.variant.subset.vcf.gz" \
  -Oz -o "${OUTDIR}/${SCAF}.variant.subset.miss80.vcf.gz"

crun.bcftools bcftools view --threads 10 \
  -i 'F_MISSING<=0.2' \
  "${OUTDIR}/${SCAF}.invariant.subset.vcf.gz" \
  -Oz -o "${OUTDIR}/${SCAF}.invariant.subset.miss80.vcf.gz"

crun.bcftools bcftools index "${OUTDIR}/${SCAF}.variant.subset.miss80.vcf.gz"
crun.bcftools bcftools index "${OUTDIR}/${SCAF}.invariant.subset.miss80.vcf.gz"

# 4) Recombine the filtered variant and invariant records into a final pixy-ready all-sites VCF.
crun.bcftools bcftools concat --threads 10 \
  --allow-overlaps \
  "${OUTDIR}/${SCAF}.variant.subset.miss80.vcf.gz" \
  "${OUTDIR}/${SCAF}.invariant.subset.miss80.vcf.gz" \
  -Ou | \
crun.bcftools bcftools sort \
  -Oz -o "${OUTDIR}/${SCAF}.pixy_ready.vcf.gz"

crun.bcftools bcftools index "${OUTDIR}/${SCAF}.pixy_ready.vcf.gz"

echo "Finished scaffold: $SCAF"

```
Runtime: 1 hr 39 min

Notes:
<br>
- `bcftools view -v snps -m2 -M2` limits variants to biallelic SNPs.
- `bcftools view -i 'ALT="."'` identifies and keeps invariants by selecting sites where the ALT field is ".".
- After sample subsetting, some formerly variant sites may lose all observed ALT alleles in the retained samples, making them invariant in the final dataset. However, these sites will still retain the ALT record from the original sample set and therefore be classified by pixy as variants. The `-a` option in `bcftools view` removes ALT alleles not seen in the retained genotypes. The record is kept but ALT is set to ".", properly reclassifying the as invariant.
- Also, the `bcftools view` command is the exception within the bcftools suite that **does** update INFO/AC and INFO/AN after subsetting with -S (unless you explicitly use the `-I/--no-update` flag), so the AC and AN info fields are updated as well in the final sample set.
- The `--recode-INFO-all` flag in the `vcftools` commands ensures that INFO fields are kept in the output vcf (these are otherwise removed by `--recode` by default).

<br>

Run some QC checks on one chromosome:
```bash
SCAF=NC_089320.1_Pverrucosa
OUTDIR=/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf_allsites/filtered

# number of records in final VCF
crun.bcftools bcftools index -n "${OUTDIR}/${SCAF}.pixy_ready.vcf.gz"

# confirm mix of invariant and variant sites
crun.bcftools bcftools query -f '%ALT\n' "${OUTDIR}/${SCAF}.pixy_ready.vcf.gz" | \
sort | uniq -c | head

# check for duplicate positions
crun.bcftools bcftools query -f '%CHROM\t%POS\n' "${OUTDIR}/${SCAF}.pixy_ready.vcf.gz" | \
sort | uniq -d | head

# count FILTER states in final file
crun.bcftools bcftools query -f '%FILTER\n' "${OUTDIR}/${SCAF}.pixy_ready.vcf.gz" | \
sort | uniq -c | sort -nr
```
```
# number of records
14488018

# breakdown of sites by variant type
14315858 .
  49683 A
  36316 C
  36331 G
  49830 T

# no duplicate positions

# breakdown of sites by FILTER state
14335363 PASS
  152655 LowQual
```
<br>

Check how many nonvariant sites were removed by depth filters:
```bash
PREDIR=/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf_allsites/

# variant sites
crun.bcftools bcftools view -H -v snps -m2 -M2 "${PREDIR}/${SCAF}_allsites_genotypes.vcf.gz" | wc -l; \
crun.bcftools bcftools view -v snps -m2 -M2 "${PREDIR}/${SCAF}_allsites_genotypes.vcf.gz" -Ou | crun.bcftools bcftools filter -e 'INFO/DP < 792 || INFO/DP > 25286' -Ou | crun.bcftools bcftools view -H | wc -l

# invariant sites
crun.bcftools bcftools view -H -i 'ALT="."' "${PREDIR}/${SCAF}_allsites_genotypes.vcf.gz" | wc -l; \
crun.bcftools bcftools view -i 'ALT="."' "${PREDIR}/${SCAF}_allsites_genotypes.vcf.gz" -Ou | crun.bcftools bcftools filter -e 'INFO/DP < 792 || INFO/DP > 25286' -Ou | crun.bcftools bcftools view -H | wc -l
```
```
1700484 # variant sites pre-depth-filtering
1560077 # variant sites post-depth-filtering

20267886 # invariant sites pre-depth-filtering
16842347 # invariant sites post-depth-filtering
```
Summary:
- Variant SNPs: 8.3% removed by DP filtering
- Invariant sites: 16.9% removed by DP filtering

So depth filtering removed more invariant sites on this chromosome than variant sites, but not drastically so. Still lots of invariant sites remaining.

<br>
<br>

## Run pixy

We'll run three separate pixy commands that address three different groups of biological questions:
- Whole dataset within-population only (summary of diversity in the species at the region scale)
- Per-location within-location summary + pairwise-between-location comparisons
- Per-island within-island summary + pairwise-between-island comparison

<br>

`pixy_10kb_array.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=pixy_10kb_array_2026-03-26
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --array=1-52%24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=5-00:00:00

set -euo pipefail

module load pixy/2.0.0.beta14

BASEDIR=/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all
VCFDIR=$BASEDIR/vcf_allsites/filtered
POPDIR=$BASEDIR/vcf
OUTBASE=$BASEDIR/pixy
SCAFLIST=/cm/shared/courses/dbarshis/barshislab/jtoy/references/genomes/pocillopora_verrucosa/ncbi_dataset/data/GCF_036669915.1/genome_regions.list
WINDOW=10000
NCORES=${SLURM_CPUS_PER_TASK}

mkdir -p $OUTBASE/{all,location,island}

SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SCAFLIST")
VCF=${VCFDIR}/${SCAF}.pixy_ready.vcf.gz

echo "Scaffold: $SCAF"
echo "VCF: $VCF"

# 1) whole dataset within-pop stats
crun.pixy pixy \
  --stats pi watterson_theta tajima_d \
  --vcf "$VCF" \
  --populations ${POPDIR}/pixy_all.pop.tsv \
  --window_size ${WINDOW} \
  --n_cores ${NCORES} \
  --output_folder ${OUTBASE}/all \
  --output_prefix ${SCAF}.all

# 2) location-level within + pairwise between-location stats
crun.pixy pixy \
  --stats pi watterson_theta tajima_d fst dxy \
  --vcf "$VCF" \
  --populations ${POPDIR}/pixy_location.pop.tsv \
  --window_size ${WINDOW} \
  --n_cores ${NCORES} \
  --fst_type hudson \
  --output_folder ${OUTBASE}/location \
  --output_prefix ${SCAF}.location

# 3) island-level within + pairwise between-island stats
crun.pixy pixy \
  --stats pi watterson_theta tajima_d fst dxy \
  --vcf "$VCF" \
  --populations ${POPDIR}/pixy_island.pop.tsv \
  --window_size ${WINDOW} \
  --n_cores ${NCORES} \
  --fst_type hudson \
  --output_folder ${OUTBASE}/island \
  --output_prefix ${SCAF}.island
```
Runtime : 26 minutes

<br>

This job returned a "Failed, Mixed, ExitCode [0-1]" status. Upon inspecting the .err files, this is because some scaffold jobs had "no invariant sites (ALT = ".")". These scaffolds were all non-chromosome scaffolds. All 14 chromosome jobs completed successfully, along with 13 non-chromosome scaffold jobs.

In addition, 9 of the non-chromosome scaffolds that completed successfully did not output a `_fst.txt` file. This is probably because these scaffolds did not have any usable **variant** sites with which Fst is normally calculated. To confirm this, run the following script to count SNPs in each scaffold that does not have a `_fst.txt` file:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/pixy/location

echo -e "SCAF\tSNPs"

for f in *_pi.txt; do
    base=${f%_pi.txt}
    
    # skip if fst file exists
    [[ -f "${base}_fst.txt" ]] && continue
    
    scaf=${base%.location}
    VCF=/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf_allsites/filtered/${scaf}.pixy_ready.vcf.gz
    
    nsnps=$(crun.bcftools bcftools view -H -v snps "$VCF" | wc -l)
    
    echo -e "${scaf}\t${nsnps}"
done | column -t
```
```
SCAF    SNPs
NW_027078169.1_Pverrucosa  0
NW_027078170.1_Pverrucosa  0
NW_027078172.1_Pverrucosa  0
NW_027078173.1_Pverrucosa  0
NW_027078175.1_Pverrucosa  0
NW_027078176.1_Pverrucosa  0
NW_027078177.1_Pverrucosa  0
NW_027078184.1_Pverrucosa  0
NW_027078186.1_Pverrucosa  0
```


This script can also be modified to summarize all scaffolds:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/pixy/location

echo -e "SCAF\tFST_file\tSNPs\tInvariants"

for f in *_pi.txt; do
    base=${f%_pi.txt}
    scaf=${base%.location}
    
    VCF=/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf_allsites/filtered/${scaf}.pixy_ready.vcf.gz
    
    if [[ -f "${base}_fst.txt" ]]; then
        fst="YES"
    else
        fst="NO"
    fi
    
    nsnps=$(crun.bcftools bcftools view -H -v snps "$VCF" | wc -l)
    ninv=$(crun.bcftools bcftools view -H -i 'ALT="."' "$VCF" | wc -l)
    
    echo -e "${scaf}\t${fst}\t${nsnps}\t${ninv}"
done | column -t
```
```
SCAF    FST_file        SNPs    Invariants
NC_089312.1_Pverrucosa     YES  271452  20172532
NC_089313.1_Pverrucosa     YES  141161  12477459
NC_089314.1_Pverrucosa     YES  241576  18034983
NC_089315.1_Pverrucosa     YES  217468  18140369
NC_089316.1_Pverrucosa     YES  117313  9576545
NC_089317.1_Pverrucosa     YES  138045  10829275
NC_089318.1_Pverrucosa     YES  126604  13058179
NC_089319.1_Pverrucosa     YES  172642  13694399
NC_089320.1_Pverrucosa     YES  172160  14315858
NC_089321.1_Pverrucosa     YES  149826  11991088
NC_089322.1_Pverrucosa     YES  93692   8959628
NC_089323.1_Pverrucosa     YES  111098  10484937
NC_089324.1_Pverrucosa     YES  110014  9872963
NC_089325.1_Pverrucosa     YES  151880  12937449
NW_027078166.1_Pverrucosa  YES  6       4663
NW_027078169.1_Pverrucosa  NO   0       335
NW_027078170.1_Pverrucosa  NO   0       1016
NW_027078171.1_Pverrucosa  YES  3       3009
NW_027078172.1_Pverrucosa  NO   0       416
NW_027078173.1_Pverrucosa  NO   0       546
NW_027078175.1_Pverrucosa  NO   0       1322
NW_027078176.1_Pverrucosa  NO   0       2643
NW_027078177.1_Pverrucosa  NO   0       564
NW_027078184.1_Pverrucosa  NO   0       16
NW_027078186.1_Pverrucosa  NO   0       263
NW_027078193.1_Pverrucosa  YES  7       1270
NW_027078194.1_Pverrucosa  YES  1       3618
```
So yes, the 9 scaffolds without a `_fst.txt` output file all had 0 variant sites.


## Filter, summarize, and plot pixy results

Import pixy output files into R and analyze each pixy run and each stat individually.
`plot_pixy_results.R`
```

```

Distribution of callable sites per 10kb window (uniform y-axis scale):
![alt text](image-3.png)

Same plot but allowing free y-axis:
![alt text](image-2.png)

Distribution of pi across all retained 10kb windows:
![alt text](image-4.png)

Faceted by chromosome:
![alt text](image-5.png)

Pi of 10kb windows across the genome:
![alt text](image-6.png)
#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name=variants_recalibrate
#SBATCH --output=/home/user/project/logs/variants_recalibrate/variants_recalibrate_slurm.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user@email.com

module load gatk-4.2.5.0
export LANG=C
export OMP_NUM_THREADS=20
conda activate bioinfo

echo Start time: $(date)

echo gatk VariantFiltration - $(date)

gatk --java-options "-Xmx60g" VariantFiltration \
    -V /home/user/project/data/variants_join/variants_join.vcf.gz \
    --filter-expression "ExcessHet > 54.69" \
    --filter-name ExcessHet \
    -O  /home/user/project/data/variants_recalibrate/variants_join_HetFiltered.vcf.gz \
    &> /home/user/project/logs/variants_recalibrate/VariantFiltration_ExcessHet.log

echo gatk MakeSitesOnlyVcf - $(date)

gatk MakeSitesOnlyVcf \
    -I /home/user/project/data/variants_recalibrate/variants_join_HetFiltered.vcf.gz \
    -O /home/user/project/data/variants_recalibrate/genotype_sites.vcf.gz \
    &> /home/user/project/logs/variants_recalibrate/makeSitesOnlyVcf.log

echo gatk VariantRecalibrator_indels - $(date)

gatk --java-options "-Xmx60g" VariantRecalibrator \
    -V /home/user/project/data/variants_recalibrate/genotype_sites.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 92.5 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -R /home/user/project/settings/genome/Homo_sapiens_assembly38.fasta \
    -mode INDEL \
    --max-gaussians 4 \
    --resource:mills,known=false,training=true,truth=true,prior=12 /home/user/project/settings/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:axiomPoly,known=false,training=true,truth=false,prior=10 /home/user/project/settings/reference/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/user/project/settings/reference/00-All.vcf.gz \
    -O /home/user/project/logs/variants_recalibrate/indels.recal \
    --tranches-file /home/user/project/logs/variants_recalibrate/indels.tranches \
    --rscript-file  /home/user/project/logs/variants_recalibrate/indels_plots.r \
    &> /home/user/project/logs/variants_recalibrate/indels_VariantRecalibrator.log

echo gatk VariantRecalibrator_snps - $(date)

gatk --java-options "-Xmx60g" VariantRecalibrator  \
    -V /home/user/project/data/variants_recalibrate/genotype_sites.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 92.5 -tranche 90.0 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
     -R /home/user/project/settings/genome/Homo_sapiens_assembly38.fasta \
    -mode SNP \
    --max-gaussians 6 \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/user/project/settings/reference/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 /home/user/project/settings/reference/1000G_omni2.5.hg38.vcf.gz  \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/user/project/settings/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 /home/user/project/settings/reference/00-All.vcf.gz \
    -O /home/user/project/logs/variants_recalibrate/snps.recal \
    --tranches-file /home/user/project/logs/variants_recalibrate/snps.tranches \
    --rscript-file  /home/user/project/logs/variants_recalibrate/snps_plots.r \
     &> /home/user/project/logs/variants_recalibrate/snps_VariantRecalibrator.log

echo gatk ApplyVQSR_indels - $(date)

gatk --java-options "-Xmx60g" \
    ApplyVQSR \
    -R /home/user/project/settings/genome/Homo_sapiens_assembly38.fasta \
    -V /home/user/project/data/variants_recalibrate/variants_join_HetFiltered.vcf.gz \
    -L /home/user/project/settings/coordinates/scattered-intervals-hg38_unix.list \
    --recal-file /home/user/project/logs/variants_recalibrate/indels.recal \
    --tranches-file //home/user/project/logs/variants_recalibrate/indels.tranches \
    --truth-sensitivity-filter-level 95.0 \
    --create-output-variant-index true \
    -mode INDEL \
    -O /home/user/project/logs/variants_recalibrate/genotypePartialRecalibrated.vcf.gz \
    &> /home/user/project/logs/variants_recalibrate/indels_ApplyVQSR.log

echo gatk ApplyVQSR_snps - $(date)

gatk --java-options "-Xmx60g" \
    ApplyVQSR \
    -R /home/user/project/settings/genome/Homo_sapiens_assembly38.fasta \
    -V /home/user/project/logs/variants_recalibrate/genotypePartialRecalibrated.vcf.gz \
    -L /home/user/project/settings/coordinates/scattered-intervals-hg38_unix.list \
    --recal-file /home/user/project/logs/variants_recalibrate/snps.recal \
    --tranches-file /home/user/project/logs/variants_recalibrate/snps.tranches \
    --truth-sensitivity-filter-level 97.5 \
    --create-output-variant-index true \
    -mode SNP \
    -O /home/user/project/data/variants_recalibrate/exomeGenotypes.vcf.gz \
    &> /home/user/project/logs/variants_recalibrate/snps_ApplyVQSR.log

echo gatk VariantsToTable - $(date)

gatk VariantsToTable \
    -V /home/user/project/data/variants_recalibrate/exomeGenotypes.vcf \
    -F CHROM -F POS -F ID -F REF -F ALT -GF GT \
    -O /home/user/project/data/variants_recalibrate/exomeGenotypes_VariantsToTable.tsv


echo End time: $(date)

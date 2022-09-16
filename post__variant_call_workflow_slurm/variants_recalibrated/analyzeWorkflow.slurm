#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name=analyzeWorkflow
#SBATCH --output=/home/user/project/logs/variants_recalibrate/analyzeWorkflow_slurm.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user@email.com

module load gatk-4.2.5.0
module load bcftools-1.15
module load htslib-1.15
export LANG=C
export OMP_NUM_THREADS=20
conda activate bioinfo

echo Start: $(date)

outputDir="/home/user/project/data/variants_recalibrate/"
logDir="/home/user/project/logs/variants_recalibrate/"
bedPath="/home/user/project/settings/coordinates/exons_ensembl106.bed"
exomeVCF="/home/user/project/data/variants_recalibrate/exomeGenotypes.vcf.gz"

echo Create BED files - $(date)

for CHROM in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
	cat $bedPath | grep -P "$CHROM\t" > $logDir/exome_$CHROM.bed
done

echo Separate VCF per CHR - $(date)


parallel="parallel --delay .2 -j 10 --joblog $logDir/analyzeWorkflow_separateChr_parallel.log "
$parallel /opt/bcftools-1.10.2/bcftools view \
	$exomeVCF \
	-R $logDir/exome_{}.bed \
	-o $outputDir/exome_{}.vcf.gz \
	-O z \
	&> $logDir/analyzeWorkflow_bcfView.log \
	 ::: chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
	chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

echo Index VCFs - $(date)

parallel tabix \
	$outputDir/exome_{}.vcf.gz ::: chr1 \
	chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

echo CollectVariantCallingMetrics per chromosome - $(date)

export OMP_NUM_THREADS=6
parallel="parallel --delay .2 -j 6 --joblog $logDir/analyzeWorkflow_chrCVCM_parallel.log "
$parallel gatk --java-options "-Xmx10g" \
	CollectVariantCallingMetrics \
	-I $outputDir/exome_{}.vcf.gz \
	--DBSNP /home/user/project/settings/reference/00-All.vcf.gz \
	-SD /home/user/project/settings/genome/Homo_sapiens_assembly38.dict \
	-O $outputDir/collectVariantCallingMetrics_{}_report \
	&> $logDir/CollectVariantCallingMetrics_chr.log \
	:::  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

echo gatk CollectVariantCallingMetrics all - $(date)

gatk CollectVariantCallingMetrics \
    -I $exomeVCF \
    --DBSNP /home/user/project/settings/reference/00-All.vcf.gz \
    -SD /home/user/project/settings/genome/Homo_sapiens_assembly38.dict \
    -O $outputDir/report_CollectVariantCallingMetrics__all_report \
    &> $logDir/CollectVariantCallingMetrics_all.log


echo Clean - $(date)

rm $logDir/exome_*.bed
rm $outputDir/exome_*.vcf.gz

echo End: $(date)

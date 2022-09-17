# Supplemental Material repository for the thesis *Caracterización de variantes genéticas y frecuencias alélicas vinculadas a hiperlipidemias a partir de genomas costarricenses* and iterations of it

This repository contains the scripts developed and used for the thesis labeled *Caracterización de variantes genéticas y frecuencias alélicas vinculadas a hiperlipidemias a partir de genomas costarricenses* and iterations of it. 

The scripts stored in the *variant_call_workflow* were developed to simplify the process of writing scripts to process genomic data in through the different stages of the bioinformatics workflow described in the thesis. The scripts ensure that log files and genomic data are is generated and documented in an orderly fashion through the workflow.

These scripts were used to be executed in the CENAT - Kabre and/or the UCR CICIMA computational cluster.

## **/thesis__variant_call_workflow_pbs**
The *thesis_variant_call_workflow_pbs* contains a set of scripts used for the thesis. In this case, scripts generates PBS scripts. It contains the following scripts:
### reads_extract.py

This script creates .sh, job lists and .pbs files to run in parallel the extraction of genomic data from CRAM files found in a specific set of genomic coordinates (should be provided) and converts it to BAM format.

*How to run from bash command line?*

> python reads_extract.py \<folder where to save scripts\> \<folder with the CRAM files to be processed\> \<3 column .bed with coordinates of interest\> \<folder where to store log files\> \<folder where to store output .bam\> \<path to the sam tools module\>

 ### 1- **reads_recalibrate.py**

This script creates .sh, job lists and .pbs files to run in parallel MarkDuplicates, BaseRecalibrator,ApplyBQSR and AnalyzeCovariates functions from GATK to the .bam files generated in last stage.

*How to run from bash command line?*

> python reads_recalibrate.py \<folder where to save scripts\> \<folder where un-recalibrated .bam files are\> \<folder where to store log files\> \<folder where to save recalibrated .bam files\> \<.fasta file corresponding to the reference genome\> \<.vcf.gz file from dbSNP\> \<path to GATK module\>

 ### 2 - **variants_individual.py**

This script creates .sh, job lists and .pbs files to run in parallel Haplotypecaller function from GATK to the recalibrated .bam files generated in last stage.

*How to run from bash command line?*

> python variants_individual.py \<folder where to save scripts\> \<folder where the recalibrated .bam files are\> \<folder where to store log files\> \<folder where save .gvcf\>  \<.fasta file corresponding to the reference genome\> \<.list with coordinates of interest\> \<path to GATK module\> \<path to java module\>

*Note: the list files should have the following format per line.*

> chr12:2000-340000

### 3 - **variants_joint.py**

This script creates .sh, job lists and .pbs files to run in parallel GenotypeGVCFs, GenomicsDBImport and MergeVcfs functions from GATK to the .gvcf files generated in last stage.

*How to run from bash command line?*

> python variants_joint.py \<folder where to save scripts\> \<folder where the .gvcf files are\> \<folder where to store log files\> \<folder where save the merged .gvcf file\> \<.fasta file corresponding to the reference genome\> \<.list with coordinates of interest\> \<.vcf.gz file from dbSNP\> \<path to GATK module\>

### 4 - **recalibrateVariants.pbs**

This is a .pbs script. It does not run in parallel the recalibration of the merged .gvcf file. In contrast to the .py scripts, parameters can't be specified to it through the command line.

### 5 - **analyzeWorkflow.pbs**

This is a .pbs script. It generates metrics that describe the output of the variant call workflow. In contrast to the .py scripts, parameters can't be specified to it through the command line.

## **/post__variant_call_workflow_pbs**
Its an iteration of the scripts found in *thesis_variant_call_workflow_pbs*.  In this case, scripts generates SLURM scripts. It contains the following scripts:

### 1 - **reads_extract.py**

This script creates .sh, job lists and .pbs files to run in parallel the extraction of genomic data from CRAM files found in a specific set of genomic coordinates (should be provided) and converts it to BAM format.

*How to run from bash command line?*

> python reads_extract.py \<folder where to save scripts\> \<.txt with list of CRAM files to be analyzed\> \<3 column .bed with coordinates of interest\> \<folder where to store log files\> \<folder where to store output .bam\> \<path to the sam tools module\> \<the amount of batches in which the samples will be processed\>

 ### 2 - **reads_recalibrate.py**

This script creates .sh, job lists and .pbs files to run in parallel MarkDuplicates, BaseRecalibrator,ApplyBQSR and AnalyzeCovariates functions from GATK to the .bam files generated in last stage.

*How to run from bash command line?*

> python reads_recalibrate.py \<folder where to save scripts\> \<.txt with list of un-recalibrated BAM files to be analyzed\> \<folder where to store log files\> \<folder where to save recalibrated .bam files\> \<.fasta file corresponding to the reference genome\> \<.vcf.gz file from dbSNP\> \<the amount of batches in which the samples will be processed\>

 ### 3 - **variants_individual.py**

This script creates .sh, job lists and .pbs files to run in parallel Haplotypecaller function from GATK to the recalibrated .bam files generated in last stage.

*How to run from bash command line?*

> python variants_individual.py \<folder where to save scripts\> <.txt with list of recalibrated BAM files to be analyzed\> \<folder where to store log files\> \<folder where save .gvcf\>  \<.fasta file corresponding to the reference genome\> \<the amount of batches in which the samples will be processed\>

*Note: the list files should have the following format per line.*

> chr12:2000-340000

### 4 - **variants_joint.py**

This script creates .sh, job lists and .pbs files to run in parallel GenotypeGVCFs, GenomicsDBImport and MergeVcfs functions from GATK to the .gvcf files generated in last stage.

*How to run from bash command line?*

> python variants_joint.py \<folder where to save scripts\>  <.txt with list of .gvcf files to be analyzed\> \<folder where to store log files\> \<folder where save the merged .gvcf file\> \<.fasta file corresponding to the reference genome\> \<.list with coordinates of interest\> \<.vcf.gz file from dbSNP\> \<the amount of batches in which the samples will be processed\>

### 5 - **recalibrateVariants.slurm**

This is a .slurm script. It does not run in parallel the recalibration of the merged .gvcf file. In contrast to the .py scripts, parameters can't be specified to it through the command line. Because of this, the script needs to be edited manually in case the script is repurposed to analyze the data found in a different environment.

___
Trabajo Final de Graduación Licenciatura en Biología Molecular y Biotecnología de Juan Carlos Valverde-Hernández, 2020-2022
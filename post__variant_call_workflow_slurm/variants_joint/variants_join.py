#!/usr/bin/python


### IMPORT

import os
import os.path
import sys
import random

### GLOBAL VARIABLES

# Location of the path that will contain the scripts for this phase
# i.e. /home/user/projects/scripts/variants_join/
scriptsDir = sys.argv[1]
# Location of folder containing the GVCF files to be genotyped and unified
# i.e. /home/user/projects/data/variants_individual/
gvcfListPath = sys.argv[2]
# Path to store log and intermediate files
# i.e. /home/user/projects/logs/variants_join/
logPath = sys.argv[3]
# Storage location of the genotyped files
# i.e. /home/user/projects/data/variants_join/
dataSavePath = sys.argv[4]
# Path to the reference genome in format. In the same directory its dictionary
# and index (fai) need to be stored
# i.e. /home/user/projects/study-specs/referenceGenome/Homo_sapiens_assembly38.fasta
refPath = sys.argv[5]
# Path to genomic coordinates list
# i.e. /home/user/projects/study-specs/genomicCoordinates/wg.list
coordinatesListPath = sys.argv[6]
# i.e. /home/user/projects/study-specs/genomicCoordinates/dbsnp.vcf
variantsDatabasePath = sys.argv[7]
# i.e. coordinate files
batchAmount = sys.argv[8]

# Location of the GVCF file map
# i.e. /home/user/projects/logs/variants_join/gvcfMap.tsv
gvcfMapPath = logPath + "/gvcfMap.tsv"
# Location of the folder containing the genomincs database import temporal directory
gdbi_dir_temp = dataSavePath + "/gdbi_temp/"
# Location of the folder containing the genomincs database import workspace directory
gdbi_dir_wp = dataSavePath + "/gdbi_workspace/"
# Location of the folder containing the genomincs database import workspace directory
ggvcf_dir = dataSavePath + "/ggvcf_partial_vcf/"
# Location of the folder containing the genomincs database import workspace directory
mgvcf2merge = logPath + "/ggvcf_partials.list"

##MAKE DIRECTORIES
if (not os.path.isdir(gdbi_dir_temp)):
    print(gdbi_dir_temp + " being created")
    os.mkdir(gdbi_dir_temp)

if (not os.path.isdir(gdbi_dir_wp)):
    print(gdbi_dir_wp + " being created")   
    os.mkdir(gdbi_dir_wp)
if (not os.path.isdir(ggvcf_dir)):
    print(ggvcf_dir + " being created")   
    os.mkdir(ggvcf_dir)


##IDENTIFY SAMPLES
gvcfListTemp = open(gvcfListPath, "r").read()
gvcfFilesPath = gvcfListTemp.split("\n")
sampleInfo = []
for i in gvcfFilesPath:
    if i == "":
        pass
    else:
        info = []
        fileName = i.split('/')
        info.append(i)
        sampleName = fileName[-1].split('.')
        sampleName = sampleName[0].split('_')
        info.append(sampleName[0])
        info.append(i)
        sampleInfo.append(info)

##IDENTIFY COORDINATES
coordinatesListTemp = open(coordinatesListPath, "r").read()
coordinatesList = coordinatesListTemp.split("\n")
random.shuffle(coordinatesList)

### FUNCTION DEFINITIONS

# Generate a map of the samples
def mapGVCF(sampleInfo):
    gvcfMapTsv = open(gvcfMapPath, "w")
    for i in sampleInfo:
        gvcfMapTsv.write(i[1] + "\t" + i[0] + "\n")
    gvcfMapTsv.close()

#Partition coordinates list
def separateCoordinates(coordinatesList, batchAmount):
    sliceLen = len(coordinatesList) // int(batchAmount)
    residue = len(coordinatesList) % int(batchAmount)
    for i in range(1, int(batchAmount) + 1 + residue ** 0):
        batch = coordinatesList[((i - 1) * sliceLen):(i * sliceLen)]
        smallCoordinatesPath = logPath + "/coordinates_" + str(i) + ".list"
        smallCoordinatesFile = open(smallCoordinatesPath, "w")
        for e in batch:
            smallCoordinatesFile.write(e + "\n")
        smallCoordinatesFile.close()


# Generate a slurm file for the scripts mentioned in the JobListFile
def gdbiSlurmFile():
    residue = len(coordinatesList) % int(batchAmount)
    for i in range(1, int(batchAmount) + 1 + residue ** 0):
        smallCoordinatesPath = logPath + "/coordinates_" + str(i) + ".list"

        slurmFile = open(scriptsDir + "/variants_join_gdbi_" + str(i) + ".slurm", "w")

        slurmFile.write("#!/bin/bash\n\n")
        slurmFile.write("#SBATCH --partition=parallel\n")
        slurmFile.write("#SBATCH --time=100:00:00\n")
        slurmFile.write("#SBATCH --nodes=1\n")    
        slurmFile.write("#SBATCH --ntasks-per-node=20\n")
        slurmFile.write("#SBATCH --job-name=variantsJoin_gdbi_" + str(i) + "\n")
        slurmFile.write("#SBATCH --output=" + logPath + "/variantsJoin_gdbi_slurm_" + str(i) + ".log\n")
        slurmFile.write("#SBATCH --mail-type=END,FAIL\n")
        slurmFile.write("#SBATCH --mail-user=user@email.com\n\n")

        slurmFile.write("module load gatk-4.2.5.0 \n")
        slurmFile.write("export OMP_NUM_THREADS=20\n")
        slurmFile.write("export LANG=C \n")
        slurmFile.write("conda activate bioinfo\n\n")

        slurmFile.write("echo Start time: $(date)\n")
        slurmFile.write("cat " + smallCoordinatesPath + " | parallel mkdir " + gdbi_dir_temp + "/gdbi_temp_{} \n")
        slurmFile.write("parallel=\"parallel --delay .1 -j 5 --joblog " + logPath + "/variants_join_gdbi_parallel_" + str(i) + ".log \"\n")
        slurmFile.write("cat " + smallCoordinatesPath + " | $parallel gatk --java-options \"-Xmx11g\" GenomicsDBImport --reader-threads 4 --genomicsdb-shared-posixfs-optimizations true --batch-size 75 --overwrite-existing-genomicsdb-workspace --sample-name-map " + gvcfMapPath + " -L {} --genomicsdb-workspace-path " + gdbi_dir_wp + "/gdbi_workspace_{}/ --tmp-dir " + gdbi_dir_temp + "/gdbi_temp_{}/ &> " + logPath + "/variants_join_gdbi_" + str(i) + ".log \n")
        slurmFile.write("echo End time: $(date)\n")

        slurmFile.close()

def ggvcfSlurmFile():
    residue = len(coordinatesList) % int(batchAmount)
    for i in range(1, int(batchAmount) + 1 + residue ** 0):
        smallCoordinatesPath = logPath + "/coordinates_" + str(i) + ".list"

        slurmFile = open(scriptsDir + "/variants_join_ggvcf_" + str(i) + ".slurm", "w")
        slurmFile.write("#!/bin/bash\n\n")
        slurmFile.write("#SBATCH --partition=parallel\n")
        slurmFile.write("#SBATCH --time=100:00:00\n")
        slurmFile.write("#SBATCH --nodes=1\n")    
        slurmFile.write("#SBATCH --ntasks-per-node=20\n")
        slurmFile.write("#SBATCH --job-name=variants_join_ggvcf_" + str(i) + "\n")
        slurmFile.write("#SBATCH --output=" + logPath + "/variants_join_ggvcf_slurm_" + str(i) + ".log\n")
        slurmFile.write("#SBATCH --mail-type=END,FAIL\n")
        slurmFile.write("#SBATCH --mail-user=user@email.com\n\n")

        slurmFile.write("module load gatk-4.2.5.0 \n")
        slurmFile.write("export OMP_NUM_THREADS=20\n")
        slurmFile.write("export LANG=C \n\n")
        slurmFile.write("conda activate bioinfo\n\n")

        slurmFile.write("echo Start time: $(date)\n")
        slurmFile.write("parallel=\"parallel --delay .2 -j 10 --joblog " + logPath + "/variants_join_ggvcf_parallel_" + str(i) + ".log \"\n")
        slurmFile.write("cat " + smallCoordinatesPath + " | $parallel gatk --java-options \"-Xmx11g\" GenotypeGVCFs -L {} -R " + refPath + " -V gendb://" + gdbi_dir_wp + "/gdbi_workspace_{}/ -D " + variantsDatabasePath + " -O " + ggvcf_dir + "/ggvcf_{}.vcf.gz &> " + logPath + "variants_join_ggvcf_" + str(i) + ".log \n")
        slurmFile.write("echo End time: $(date)\n")

    slurmFile.close()

def mgcvfSlurmFile():
    slurmFile = open(scriptsDir + "/variants_join_mgcvf.slurm", "w")
    slurmFile.write("#!/bin/bash\n\n")
    slurmFile.write("#SBATCH --partition=parallel\n")
    slurmFile.write("#SBATCH --time=100:00:00\n")
    slurmFile.write("#SBATCH --nodes=1\n")    
    slurmFile.write("#SBATCH --ntasks-per-node=20\n")
    slurmFile.write("#SBATCH --job-name=variants_join_mgcvf\n")
    slurmFile.write("#SBATCH --output=" + logPath + "/variants_join_mgcvf_slurm.log\n")
    slurmFile.write("#SBATCH --mail-type=END,FAIL\n")
    slurmFile.write("#SBATCH --mail-user=user@email.com\n\n")

    slurmFile.write("module load gatk-4.2.5.0 \n")
    slurmFile.write("export LANG=C \n\n")
    slurmFile.write("conda activate bioinfo\n\n")

    slurmFile.write("echo Start time: $(date)\n")
    slurmFile.write("ls " + ggvcf_dir + "/ggvcf_*.vcf.gz > " + mgvcf2merge + "\n")
    slurmFile.write("gatk MergeVcfs -I " + mgvcf2merge + " -O " + dataSavePath + "/variants_join.vcf.gz &> " + logPath + "/variants_join_mgcvf.log\n")
    slurmFile.write("echo End time: $(date)\n")

    slurmFile.close()

### IMPLEMENTATION

# Create slurm files

mapGVCF(sampleInfo)

separateCoordinates(coordinatesList, batchAmount)

gdbiSlurmFile()

ggvcfSlurmFile()

mgcvfSlurmFile()

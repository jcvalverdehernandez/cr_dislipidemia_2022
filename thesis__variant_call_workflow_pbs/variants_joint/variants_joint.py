#!/usr/bin/python


### IMPORT

import os
import os.path
import sys

### GLOBAL VARIABLES

# Location of the path that will contain the scripts for this phase
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/
scriptsDir = sys.argv[1]
# Location of folder containing the GVCF files to be genotyped and unified
# i.e. /home/jcvh/cr_dislipidemia/data/03_recalibrateBAM/
gvcfDir = sys.argv[2]
# Path to store log and intermediate files
# i.e. /home/jcvh/cr_dislipidemia/logs/04_createGenotype/
logPath = sys.argv[3]
# Storage location of the genotyped files
# i.e. /home/jcvh/cr_dislipidemia/data/04_createGenotype/
dataSavePath = sys.argv[4]
# Path to the reference genome in format. In the same directory its dictionary
# and index (fai) need to be stored
# i.e. /home/jcvh/cr_dislipidemia/study-specs/referenceGenome/Homo_sapiens_assembly38.fasta
refPath = sys.argv[5]
# Path to genomic coordinates list
# i.e. /home/jcvh/cr_dislipidemia/study-specs/genomicCoordinates/wg.list
coordinatesListPath = sys.argv[6]
# i.e. /home/jcvh/cr_dislipidemia/study-specs/genomicCoordinates/dbsnp.vcf
variantsDatabasePath = sys.argv[7]
# Name of the GATK module
# i.e. /opt/gatk-4.2.0.0/gatk
gatkPath = sys.argv[8]

# Identifier of the samples in the gvcfDir folder
gvcfPointersDict_Temp= os.listdir(gvcfDir)
gvcfPointersDict= {}
for i in gvcfPointersDict_Temp:
    file_name = i.split('.')
    if (file_name[-2] == "gvcf") or (file_name[-2] == "vcf"):
        gvcfPointersDict[file_name[0]] = gvcfDir + '.'.join(file_name)
    else:
        pass

# Location of the GVCF file map
# i.e. /home/jcvh/cr_dislipidemia/logs/04_createGenotype/gvcfMap.tsv
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
    os.mkdir(gdbi_dir_temp)
if (not os.path.isdir(gdbi_dir_wp)):
    os.mkdir(gdbi_dir_wp)
if (not os.path.isdir(ggvcf_dir)):
    os.mkdir(ggvcf_dir)

### FUNCTION DEFINITIONS

# Generate a map of the
def mapGVCF():
    gvcfMapTsv = open(gvcfMapPath, "w")
    for i in gvcfPointersDict:
        gvcfMapTsv.write(i + "\t" + gvcfPointersDict[i] + "\n")
    gvcfMapTsv.close()


# Generate a pbs file for the scripts mentioned in the JobListFile
def gdbiPBSFile():
    pbsFile = open(scriptsDir + "/variants_joint_gdbi.pbs", "w")
    pbsFile.write("#!/bin/bash\n\n")
    pbsFile.write("#PBS -V\n")
    pbsFile.write("#PBS -N variants_joint_gdbi\n")
    pbsFile.write("#PBS -q default\n")    
    pbsFile.write("#PBS -l walltime=200:00:00\n")
    pbsFile.write("#PBS -o " + logPath + "/variants_joint_gdbi_pbs.log\n")
    pbsFile.write("#PBS -l nodes=1:ppn=" + "20" + "\n")
    pbsFile.write("#PBS -M jcvalverdehernandez@gmail.com\n")
    pbsFile.write("#PBS -m abe\n\n")
    pbsFile.write("export OMP_NUM_THREADS=10\n")

    pbsFile.write("echo Start time: $(date)\n")

    pbsFile.write("cat " + coordinatesListPath + " | parallel mkdir " + gdbi_dir_temp + "/gdbi_temp_{} \n")

    pbsFile.write("parallel=\"parallel --delay .2 -j 10 --joblog " + logPath + "/variants_joint_gdbi_parallel.log \"\n")

    pbsFile.write("cat " + coordinatesListPath + " | $parallel /opt/gatk-4.2.0.0/gatk --java-options \"-Xmx6g\" GenomicsDBImport --sample-name-map " + gvcfMapPath + " -L {} --genomicsdb-workspace-path " + gdbi_dir_wp + "/gdbi_workspace_{}/ --tmp-dir " + gdbi_dir_temp + "/gdbi_temp_{}/ &> " + logPath + "/variants_joint_gdbi.log \n")

    pbsFile.write("echo End time: $(date)\n")

    pbsFile.close()

def ggvcfPBSFile():
    pbsFile = open(scriptsDir + "/variants_joint_ggvcf.pbs", "w")
    pbsFile.write("#!/bin/bash\n\n")
    pbsFile.write("#PBS -V\n")
    pbsFile.write("#PBS -N variants_joint_ggvcf\n")
    pbsFile.write("#PBS -q default\n")    
    pbsFile.write("#PBS -l walltime=200:00:00\n")
    pbsFile.write("#PBS -o " + logPath + "/variants_joint_ggvcf_pbs.log\n")
    pbsFile.write("#PBS -l nodes=1:ppn=" + "20" + "\n")
    pbsFile.write("#PBS -M jcvalverdehernandez@gmail.com\n")
    pbsFile.write("#PBS -m abe\n\n")
    pbsFile.write("export OMP_NUM_THREADS=10\n")

    pbsFile.write("echo Start time: $(date)\n")

    pbsFile.write("parallel=\"parallel --delay .2 -j 10 --joblog " + logPath + "/variants_joint_ggvcf_parallel.log \"\n")

    pbsFile.write("cat " + coordinatesListPath + " | $parallel /opt/gatk-4.2.0.0/gatk --java-options \"-Xmx6g\" GenotypeGVCFs --genomicsdb-shared-posixfs-optimizations true -L {} -R " + refPath + " -V gendb://" + gdbi_dir_wp + "/gdbi_workspace_{}/ -D " + variantsDatabasePath + " -O " + ggvcf_dir + "/ggvcf_{}.vcf.gz &> " + logPath + "variants_joint_ggvcf.log \n")

    pbsFile.write("echo End time: $(date)\n")

    pbsFile.close()

def mgcvfPBSFile():
    pbsFile = open(scriptsDir + "/variants_joint_mgcvf.pbs", "w")
    pbsFile.write("#!/bin/bash\n\n")
    pbsFile.write("#PBS -V\n")
    pbsFile.write("#PBS -N variants_joint_mgcvf\n")
    pbsFile.write("#PBS -q default\n")    
    pbsFile.write("#PBS -l walltime=200:00:00\n")
    pbsFile.write("#PBS -o " + logPath + "/variants_joint_mgcvf_pbs.log\n")
    pbsFile.write("#PBS -l nodes=1:ppn=" + "20" + "\n")
    pbsFile.write("#PBS -M jcvalverdehernandez@gmail.com\n")
    pbsFile.write("#PBS -m abe\n\n")

    pbsFile.write("echo Start time: $(date)\n")

    pbsFile.write("ls " + ggvcf_dir + "/ggvcf_*.vcf.gz > " + mgvcf2merge + "\n")

    pbsFile.write("/opt/gatk-4.2.0.0/gatk MergeVcfs -I " + mgvcf2merge + " -O " + dataSavePath + "/variants_joint.vcf.gz &> " + logPath + "/variants_joint_mgcvf.log\n")

    pbsFile.write("echo End time: $(date)\n")

    pbsFile.close()

### IMPLEMENTATION

# Create pbs files
mapGVCF()

gdbiPBSFile()

ggvcfPBSFile()

mgcvfPBSFile()

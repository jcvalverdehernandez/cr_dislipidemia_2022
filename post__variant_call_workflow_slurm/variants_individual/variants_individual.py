#!/usr/bin/python


### IMPORT

import os
import os.path
import sys

### GLOBAL VARIABLES

# Location of the path that will contain the scripts for this phase
# i.e. /home/user/project/scripts/variants_individual/
scriptsDir = sys.argv[1]
# Location of folder containing the BAM files to be genotyped and unified
# i.e. /home/user/project/data/reads_recalibrate/
bamListPath = sys.argv[2]
# Path to store log and intermediate files
# i.e. /home/user/project/logs/variants_individual/
logPath = sys.argv[3]
# Storage location of the genotyped files
# i.e. /home/user/project/data/variants_individual/
dataSavePath = sys.argv[4]
# Path to the reference genome in format. In the same directory its dictionary
# and index (fai) need to be stored
# i.e. /home/user/project/study-specs/referenceGenome/Homo_sapiens_assembly38.fasta
refPath = sys.argv[5]
# Number of batches
batchAmount = sys.argv[6]

# Identifier of the samples in the bamDir folder
bamListTemp = open(bamListPath, "r").read()
bamFilesPath = bamListTemp.split("\n")
sampleInfo = []
for i in bamFilesPath:
    if i == "":
        pass
    else:
        info = []
        fileName = i.split('/')
        info.append('/'.join(fileName))
        sampleName = fileName[-1].split('.')
        info.append('.'.join(sampleName[:-1]))
        sampleInfo.append(info)


# Location to save generated slurm script
# i.e. /home/user/project/scripts/variants_individual/shScripts/
shScriptDirName = scriptsDir + "/shScripts/"
# Location to save generated scripts list
# i.e. /home/user/project/scripts/variants_individual/jobLists
joblistsDirName = scriptsDir + "/jobLists/"
# Location to save generated slurm script
# i.e. /home/user/project/scripts/variants_individual/slurm/
slurmDirName = scriptsDir + "/slurm/"
# Location of the GVCF file map
# i.e. /home/user/project/logs/variants_individual/gvcfMap.tsv
gvcfMapPath = logPath + "gvcfMap.tsv"
#Path to sh job list
# i.e. /home/user/project/scripts/variants_individual/jobLists/jobList.txt
jobFilePath = joblistsDirName + "jobList.txt"

##MAKE DIRECTORIES
# i.e. /home/user/project/scripts/variants_individual/shScripts/
if (not os.path.isdir(shScriptDirName)):
    os.mkdir(shScriptDirName)
# i.e. /home/user/project/scripts/variants_individual/jobLists/
if (not os.path.isdir(joblistsDirName)):
    os.mkdir(joblistsDirName)
# i.e. /home/user/project/scripts/variants_individual/slurm/
if (not os.path.isdir(slurmDirName)):
    os.mkdir(slurmDirName)


### FUNCTION DEFINITIONS

shScriptPathList = []
def createScripts(info):
    # Create the shell script file
    bamFilePath = info[0]
    bamID = info[1]
    shScriptPath = shScriptDirName + bamID + ".sh"
    shScriptPathList.append(shScriptPath)
    individualGVCFFilePath = dataSavePath + bamID + ".gvcf.gz"
    haploCallLogPath = logPath + bamID + "_haploCall.log"
    sha256sumLogPath = logPath + "sha256sum.log"
    md5sumLogPath = logPath + "md5sum.log"

    shellFile = open(shScriptPath, "w")
    # Write to the file
    shellFile.write("#!/bin/bash\n\n")
    shellFile.write("start=`date +%s`\n")
    shellFile.write("start_moment=$(date)\n")
    shellFile.write("echo Operation gatk HaplotypeCaller starts $(date) " + bamID + "\n")
    shellFile.write('gatk --java-options "-Xmx4g" HaplotypeCaller -ERC GVCF -R ' + refPath + " -I " + bamFilePath + " -O " + individualGVCFFilePath + " -G StandardAnnotation -G StandardHCAnnotation &> " + haploCallLogPath + "\n")
    shellFile.write("echo Operation sha256sum starts $(date) " + bamID + "\n")
    shellFile.write("sha256sum " + individualGVCFFilePath + " >> " + sha256sumLogPath + "\n")
    shellFile.write("echo Operation md5sum starts $(date) " + bamID + "\n")
    shellFile.write("md5sum " + individualGVCFFilePath + " >> " + md5sumLogPath + "\n")
    shellFile.write("end=`date +%s`\n")
    shellFile.write("end_moment=$(date)\n")
    shellFile.write("echo Job for " + bamID + " started on $start_moment and ended at $end_moment - Total time: $((end-start)) \n")
    shellFile.close()
    os.chmod(shScriptPath, 0o744)

def makeJobListFiles(batchAmount):
    sliceLen = len(shScriptPathList) // int(batchAmount)
    residue = len(shScriptPathList) % int(batchAmount)
    for i in range(1, int(batchAmount) + 1 + residue ** 0):
        batch = shScriptPathList[((i - 1) * sliceLen):(i * sliceLen)]
        jobFilePath = joblistsDirName + "/JobList_" + str(i) + ".txt"
        jobListFile = open(jobFilePath, "w")
        for e in batch:
            jobListFile.write(e + "\n")
        jobListFile.close()

    
# Generate a slurm file for the scripts mentioned in the JobListFile
def makeslurmJobFiles(batchAmount):
    residue = len(shScriptPathList) % int(batchAmount)
    for i in range(1, int(batchAmount) + 1 + residue ** 0):
        jobFilePath = joblistsDirName + "/JobList_" + str(i) + ".txt"

        slurmFile = open(slurmDirName + "/variants_individual_" + str(i) + ".slurm", "w")
 
        slurmFile.write("#!/bin/bash\n\n")
        slurmFile.write("#SBATCH --partition=parallel\n")
        slurmFile.write("#SBATCH --time=80:00:00\n")
        slurmFile.write("#SBATCH --nodes=1\n")    
        slurmFile.write("#SBATCH --ntasks-per-node=20\n")
        slurmFile.write("#SBATCH --job-name=variants_individual_" + str(i) + "\n")
        slurmFile.write("#SBATCH --output=" + logPath + "/variants_individual_slurm_" + str(i) + ".log\n")
        slurmFile.write("#SBATCH --mail-type=END,FAIL\n")
        slurmFile.write("#SBATCH --mail-user=user@email.com\n\n")

        slurmFile.write("module load gatk-4.2.5.0\n")
        slurmFile.write("export OMP_NUM_THREADS=6\n")
        slurmFile.write("export LANG=C\n")
        slurmFile.write("conda activate bioinfo\n\n")

        slurmFile.write("echo Start time: $(date)\n")
        slurmFile.write("parallel=\"parallel --delay .2 -j 6 --joblog " + logPath + "/variants_individual_parallel_" + str(i) + ".log --resume\"\n")
        slurmFile.write("$parallel < " + jobFilePath + "\n")
        slurmFile.write("echo End time: $(date)\n")
        slurmFile.close()


### IMPLEMENTATION


bamListTemp = open(bamListPath, "r").read()
bamFilesPath = bamListTemp.split("\n")

#Create shell scripts
for i in sampleInfo:
    createScripts(i)

# Create job list
makeJobListFiles(batchAmount)

# Create slurm file
makeslurmJobFiles(batchAmount)

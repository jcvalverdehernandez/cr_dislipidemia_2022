#!/usr/bin/python


### IMPORT

import os
import os.path
import sys

### GLOBAL VARIABLES


# Location of the path that will contain the scripts for this phase
# i.e. /home/user/project/scripts/reads_recalibrate/
scriptsDir = sys.argv[1]
# Location of folder containing the BAM files to recalibrate - focused for ILLUMINA
# i.e. /home/user/project/data/reads_extract/
bamListPath = sys.argv[2]
# Path to store log and intermediate files
# i.e. /home/user/project/logs/reads_recalibrate/
logPath = sys.argv[3]
# Storage location of recalibrated BAM files
# i.e. /home/user/project/data/reads_recalibrate/
dataSavePath = sys.argv[4]
# Path to the reference genome in format. In the same directory its dictionary
# and index (fai) need to be stored
# i.e. /home/user/project/study-specs/referenceGenome/Homo_sapiens_assembly38.fasta 
refPath = sys.argv[5]
# Path to the VCF file containing the known variant sites (such as dbSNP or 1k genomes database).
# Its respective index file needs to be stored in the same folder (generated with gatk IndexFeatureFile)
# i.e. /home/user/project/study-specs/variantsReference/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
variantsPath = sys.argv[6]
# Number of batches
batchAmount = sys.argv[7]


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
# i.e. /work/user/project/scripts/reads_recalibrate/shScripts/
shScriptDirName = scriptsDir + "/shScripts/"
# Location to save generated scripts list
# i.e. /work/user/project/scripts/reads_recalibrate/jobLists
joblistsDirName = scriptsDir + "/jobLists/"
# Location to save generated slurm script
# i.e. /work/user/project/scripts/reads_recalibrate/slurm/
slurmDirName = scriptsDir + "/slurm/"


##MAKE DIRECTORIES
# i.e. /work/user/project/scripts/reads_recalibrate/shScripts/
if (not os.path.isdir(shScriptDirName)):
    os.mkdir(shScriptDirName)
# i.e. /work/user/project/scripts/reads_recalibrate/jobLists/
if (not os.path.isdir(joblistsDirName)):
    os.mkdir(joblistsDirName)
# i.e. /work/user/project/scripts/reads_recalibrate/slurm/
if (not os.path.isdir(slurmDirName)):
    os.mkdir(slurmDirName)

### FUNCTION DEFINITIONS

shScriptPathList = []
def createScripts(info):
    # Create the shell script file
    bamFilePath = info[0]
    shScriptPath = shScriptDirName + info[1] + ".sh"
    shScriptPathList.append(shScriptPath)
    bamFilePathNoDuplicates = dataSavePath + info[1] + "_nodup.bam"
    MarkDuplicatesuplicatesMetrics = logPath + info[1] + "_DuplicateMetrics.txt"
    MarkDuplicatesLogFile = logPath + info[1] + "_MarkDuplicates.log"
    baseRecalTablePath_1 = logPath + info[1] + "_baseRecal1.table"
    baseRecalLogPath_1 = logPath + info[1] + "_baseRecal1.log"
    bamRecalibratedPath = dataSavePath + info[1] + "_recal.bam"
    applyBQSRLogPath = logPath + info[1] + "_applyBQSR.log"
    baseRecalTablePath_2 = logPath + info[1] + "_baseRecal2.table"
    baseRecalLogPath_2 = logPath + info[1] + "_baseRecal2.log"
    covariantsCsvPath =logPath + info[1] + "_analyzeCovariants.csv"
    covariantsLogPath =logPath + info[1] + "_analyzeCovariants.log"
    sha256sumLogPath = logPath + "sha256sum.log"
    md5sumLogPath = logPath + "md5sum.log"
    shellFile = open(shScriptPath, "w")
    # Write to the file
    shellFile.write("#!/bin/bash\n\n")
    shellFile.write("start=`date +%s`\n")
    shellFile.write("start_moment=$(date)\n")
    shellFile.write("echo Operation gatk MarkDuplicates starts $(date) " + info[1] + "\n")
    shellFile.write("gatk MarkDuplicates -I " + bamFilePath + " -O " + bamFilePathNoDuplicates + " -M " + MarkDuplicatesuplicatesMetrics + " --REMOVE_DUPLICATES true &> " + MarkDuplicatesLogFile + "\n")
    shellFile.write("echo Operation gatk BaseRecalibrator 1 starts $(date) " + info[1] + "\n")
    shellFile.write("gatk BaseRecalibrator -I " + bamFilePathNoDuplicates + " -R " + refPath + " --known-sites " + variantsPath + " -O " + baseRecalTablePath_1 + " &> " + baseRecalLogPath_1 + "\n")
    shellFile.write("echo Operation gatk ApplyBQSR starts $(date) " + info[1] + "\n")
    shellFile.write("gatk ApplyBQSR -I " + bamFilePathNoDuplicates + " -R " + refPath + " --bqsr-recal-file " + baseRecalTablePath_1 + " -O " + bamRecalibratedPath + " &> " + applyBQSRLogPath + "\n")
    shellFile.write("echo Operation gatk BaseRecalibrator 2 starts $(date) " + info[1] + "\n")
    shellFile.write("gatk BaseRecalibrator -I " + bamRecalibratedPath + " -R " + refPath + " --known-sites " + variantsPath + " -O " + baseRecalTablePath_2 + " &> " + baseRecalLogPath_2 + "\n")
    shellFile.write("echo Operation gatk AnalyzeCovariate starts $(date) " + info[1] + "\n")
    shellFile.write("gatk AnalyzeCovariates -before " + baseRecalTablePath_1 + " -after " + baseRecalTablePath_2 + " -csv " + covariantsCsvPath + " &> " + covariantsLogPath + "\n")
    shellFile.write("echo Operation sha256sum starts $(date) " + info[1] + "\n")
    shellFile.write("sha256sum " + bamRecalibratedPath + " >> " + sha256sumLogPath + "\n")
    shellFile.write("echo Operation md5sum starts $(date) " + info[1] + "\n")
    shellFile.write("md5sum " + bamRecalibratedPath + " >> " + md5sumLogPath + "\n")
    shellFile.write("end=`date +%s`\n")
    shellFile.write("end_moment=$(date)\n")
    shellFile.write("echo Job for " + info[1] + " started on $start_moment and ended at $end_moment - Total time: $((end-start)) \n")
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

bamListTemp = open(bamListPath, "r").read()
bamFilesPath = bamListTemp.split("\n")

# Generate a slurm file for the scripts mentioned in the JobListFile
def makeslurmJobFiles(batchAmount):
    residue = len(shScriptPathList) % int(batchAmount)
    for i in range(1, int(batchAmount) + 1 + residue ** 0):
        jobFilePath = joblistsDirName + "/JobList_" + str(i) + ".txt"

        slurmFile = open(slurmDirName + "/reads_recalibrate_" + str(i) + ".slurm", "w")
        slurmFile.write("#!/bin/bash\n\n")
        slurmFile.write("#SBATCH --partition=parallel\n")
        slurmFile.write("#SBATCH --time=80:00:00\n")
        slurmFile.write("#SBATCH --nodes=1\n")    
        slurmFile.write("#SBATCH --ntasks-per-node=20\n")
        slurmFile.write("#SBATCH --job-name=\"reads_recalibrate_" + str(i) + "\"\n")
        slurmFile.write("#SBATCH --output=" + logPath + "/reads_recalibrate_slurm_" + str(i) + ".log\n")
        slurmFile.write("#SBATCH --mail-type=END,FAIL\n")
        slurmFile.write("#SBATCH --mail-user=juser@email.com\n\n")

        slurmFile.write("export LANG=C\n")
        slurmFile.write("export OMP_NUM_THREADS=4\n")
        slurmFile.write("module load gatk-4.2.5.0\n")
        slurmFile.write("module load samtools-1.15\n")
        slurmFile.write("conda activate bioinfo\n\n")

        slurmFile.write("echo Start time: $(date)\n")
        slurmFile.write("parallel=\"parallel --delay .2 -j 4 --line-buffer --joblog " + logPath + "/reads_recalibrate_parallel_" + str(i) + ".log \"\n")
        slurmFile.write("$parallel < " + jobFilePath + "\n")
        slurmFile.write("echo End time: $(date)\n")
        slurmFile.close()

### IMPLEMENTATION

#Create shell scripts
for i in sampleInfo:
    createScripts(i)

# Create job list
makeJobListFiles(batchAmount)


# Create slurm file
makeslurmJobFiles(batchAmount)

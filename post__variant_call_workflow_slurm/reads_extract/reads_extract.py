#!/usr/bin/python


### IMPORT

import os
import os.path
import sys


### GLOBAL VARIABLES

# Location of the path that will contain the scripts for this phase
# i.e. /work/user/project/scripts/02_createBAM/
scriptsDir = sys.argv[1]
# Location of folder containing the list of CRAM files
# i.e. /work/user/project/data/raw_reads/
# ! text with path of cram files to be processed
cramListPath = sys.argv[2]
# Location of the BED file that contains the regions of interest
# i.e. /work/user/project/study-specs/analysis-coordinates/roi_coordinates.bed
bedPath = sys.argv[3]
# Path to store log files
# i.e. /work/user/project/logs/02_createBAM/
logPath = sys.argv[4]
# Path for generated data files
# i.e. /work/user/project/data/02_createBAM/
dataSavePath = sys.argv[5]
# Number of batches
batchAmount = sys.argv[6]

# Identifier of the samples in the bamDir folder
cramListTemp = open(cramListPath, "r").read()
cramFilesPath = cramListTemp.split("\n")
sampleInfo = []
for i in cramFilesPath:
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
# i.e. /work/user/project/scripts/03_recalibrateBAM/shScripts/
shScriptDirName = scriptsDir + "/shScripts/"
# Location to save generated scripts
# i.e. /work/user/project/scripts/02_createBAM/jobLists
joblistsDirName = scriptsDir + "/jobLists/"
# Location to save generated slurm script
# i.e. /work/user/project/scripts/02_createBAM/slurm/
slurmDirName = scriptsDir + "/slurm/"


##MAKE DIRECTORIES
# i.e. /work/user/project/scripts/02_createBAM/shScripts
if (not os.path.isdir(shScriptDirName)):
    os.mkdir(shScriptDirName)
# i.e. /work/user/project/scripts/02_createBAM/batches/medellin
if (not os.path.isdir(joblistsDirName)):
    os.mkdir(joblistsDirName)
# i.e. /work/user/project/scripts/03_recalibrateBAM/slurm/
if (not os.path.isdir(slurmDirName)):
    os.mkdir(slurmDirName)

### FUNCTION DEFINITIONS

shScriptPathList = []

def createScripts(info):
    # Create the shell script file
    shScriptPath = shScriptDirName + info[1] + "_roi.sh"
    shScriptPathList.append(shScriptPath)
    cramFilePath = info[0]
    statsOutFile = logPath + info[1] + "_statsBAM.log"
    sha256sumLogPath = logPath + "sha256sum.log"
    md5sumLogPath = logPath + "md5sum.log"
    bamFilePath = dataSavePath + info[1] + "_roi.bam"
    bamErrPath = logPath + info[1] + "_BAM_roi.log"
    bamHeadTailSummary = logPath + info[1] + "_headTailBAM_roi.log"
    shellFile = open(shScriptPath, "w")
    # Write to the file
    shellFile.write("#!/bin/bash\n\n")
    shellFile.write("start=`date +%s`\n")
    shellFile.write("start_moment=$(date)\n")
    shellFile.write("echo Operation sha256sum starts $(date) " + info[1] + "\n")
    shellFile.write("sha256sum " + cramFilePath + " >> " + sha256sumLogPath + "\n")
    shellFile.write("echo Operation md5sum starts $(date) " + info[1] + "\n")
    shellFile.write("md5sum " + cramFilePath + " >> " + md5sumLogPath + "\n")
    shellFile.write("echo Operation samtools view starts $(date) " + info[1] + "\n")
    shellFile.write("samtools view -h -b -L " + bedPath + " " + cramFilePath + " > " + bamFilePath + " 2> " + bamErrPath + "\n")
    shellFile.write("echo Operation sha256sum view starts $(date) " + info[1] + "\n")
    shellFile.write("sha256sum " + bamFilePath + " >> " + sha256sumLogPath + "\n")
    shellFile.write("echo Operation md5sum view starts $(date) " + info[1] + "\n")
    shellFile.write("md5sum " + bamFilePath + " >> " + md5sumLogPath + "\n")
    shellFile.write("echo Operation bam head-tail starts $(date) " + info[1] + "\n")
    shellFile.write("samtools view -H " + bamFilePath + " >> " + bamHeadTailSummary + "\n")
    shellFile.write("samtools view " + bamFilePath + " | tail -n 50  >> " + bamHeadTailSummary + "\n")
    shellFile.write("echo Operation samtools stats starts $(date) " + info[1] + "\n")
    shellFile.write("samtools stats " + bamFilePath + " >> " + statsOutFile + "\n")
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

cramListTemp = open(cramListPath, "r").read()
cramFilesPath = cramListTemp.split("\n")


# Generate a SLURM file for the scripts mentioned in the JobListFile
def makeslurmJobFiles(batchAmount):
    residue = len(shScriptPathList) % int(batchAmount)
    for i in range(1, int(batchAmount) + 1 + residue ** 0):
        jobFilePath = joblistsDirName + "/JobList_" + str(i) + ".txt"

        slurmFile = open(slurmDirName + "/reads_extract_" + str(i) + ".slurm", "w")
        slurmFile.write("#!/bin/bash\n\n")
        slurmFile.write("#SBATCH --partition=parallel\n")
        slurmFile.write("#SBATCH --time=80:00:00\n")
        slurmFile.write("#SBATCH --nodes=1\n")    
        slurmFile.write("#SBATCH --ntasks-per-node=20\n")
        slurmFile.write("#SBATCH --job-name=\"reads_extract_" + str(i) + "\"\n")
        slurmFile.write("#SBATCH --output=" + logPath + "/reads_extract_slurm_" + str(i) + ".log\n")
        slurmFile.write("#SBATCH --mail-type=END,FAIL\n")
        slurmFile.write("#SBATCH --mail-user=user@email.com\n\n")
        
        slurmFile.write("module load samtools-1.15 \n")
        slurmFile.write("export LANG=C \n")
        slurmFile.write("export OMP_NUM_THREADS=10\n\n")
        slurmFile.write("parallel=\"parallel --delay .2 -j 10 --line-buffer --joblog " + logPath + "/reads_extract_parallel_" + str(i) + ".log \"\n")
        slurmFile.write("$parallel < " + jobFilePath + "\n")
        slurmFile.close()

### IMPLEMENTATION


#Create shell scripts
for i in sampleInfo:
    createScripts(i)

# Create job list
makeJobListFiles(batchAmount)

# Create slurm file
makeslurmJobFiles(batchAmount)

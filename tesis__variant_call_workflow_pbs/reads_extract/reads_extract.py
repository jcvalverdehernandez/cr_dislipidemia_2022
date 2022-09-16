#!/usr/bin/python


### IMPORT

import os
import os.path
import sys


### GLOBAL VARIABLES

# Location of the path that will contain the scripts for this phase
# i.e. /work/jvalverde/cr_dislipidemia/scripts/02_createBAM/
scriptsDir = sys.argv[1]
# Location of folder containing the list of CRAM files
# i.e. /work/jvalverde/cr_dislipidemia/data/raw_reads/
cramPath = sys.argv[2]
# Location of the BED file that contains the regions of interest
# i.e. /work/jvalverde/cr_dislipidemia/study-specs/analysis-coordinates/goi_coordinates.bed
bedPath = sys.argv[3]
# Path to store log files
# i.e. /work/jvalverde/cr_dislipidemia/logs/02_createBAM/
logPath = sys.argv[4]
# Path for generated data files
# i.e. /work/jvalverde/cr_dislipidemia/data/02_createBAM/
dataSavePath = sys.argv[5]
# Name of the batch
# Samtools module
# i.e. samtools/1.9
samtoolsModule = sys.argv[6]


# Identifier of the samples in the bamDir folder
cramDirFilesTemp = os.listdir(cramPath)
cramDirFiles = []
for i in cramDirFilesTemp:
    file_name = i.split('.')
    if file_name[-1] == "cram":
        cramDirFiles.append('.'.join(file_name[:-1]))
    else:
        pass
# Location to save generated pbs script
# i.e. /work/jvalverde/cr_dislipidemia/scripts/03_recalibrateBAM/shScripts/
shScriptDirName = scriptsDir + "/shScripts/"
# Location to save generated scripts
# i.e. /work/jvalverde/cr_dislipidemia/scripts/02_createBAM/jobLists
joblistsDirName = scriptsDir + "/jobLists/"
# Location to save generated slurm script
# i.e. /work/jvalverde/cr_dislipidemia/scripts/02_createBAM/pbs/
pbsDirName = scriptsDir + "/pbs/"


##MAKE DIRECTORIES
# i.e. /work/jvalverde/cr_dislipidemia/scripts/02_createBAM/shScripts
if (not os.path.isdir(shScriptDirName)):
    os.mkdir(shScriptDirName)
# i.e. /work/jvalverde/cr_dislipidemia/scripts/02_createBAM/batches/medellin
if (not os.path.isdir(joblistsDirName)):
    os.mkdir(joblistsDirName)
# i.e. /work/jvalverde/cr_dislipidemia/scripts/03_recalibrateBAM/pbs/
if (not os.path.isdir(pbsDirName)):
    os.mkdir(pbsDirName)

### FUNCTION DEFINITIONS

shScriptPathDict = {}

def createScripts(sampleID):
    # Create the shell script file
    shScriptPath = shScriptDirName + sampleID + "_goi.sh"
    shScriptPathDict[sampleID] = shScriptPath
    cramFilePath = cramPath + sampleID + ".cram"
    statsOutFile = logPath + sampleID + "_statsBAM.log"
    sha256sumLogPath = logPath + "sha256sum.log"
    md5sumLogPath = logPath + "md5sum.log"
    bamFilePath = dataSavePath + sampleID + "_goi.bam"
    bamErrPath = logPath + sampleID + "_BAM_goi.log"
    bamHeadTailSummary = logPath + sampleID + "_headTailBAM_goi.log"
    shellFile = open(shScriptPath, "w")
    # Write to the file
    shellFile.write("#!/bin/bash\n\n")
    shellFile.write("sha256sum " + cramFilePath + " >> " + sha256sumLogPath + "\n")
    shellFile.write("md5sum " + cramFilePath + " >> " + md5sumLogPath + "\n")
    shellFile.write(samtoolsModule + " view -h -b -L " + bedPath + " " + cramFilePath + " > " + bamFilePath + " 2> " + bamErrPath + "\n")
    shellFile.write("sha256sum " + bamFilePath + " >> " + sha256sumLogPath + "\n")
    shellFile.write("md5sum " + bamFilePath + " >> " + md5sumLogPath + "\n")
    shellFile.write(samtoolsModule + " view -H " + bamFilePath + " >> " + bamHeadTailSummary + "\n")
    shellFile.write(samtoolsModule + " view " + bamFilePath + " | tail -n 50  >> " + bamHeadTailSummary + "\n")
    shellFile.write(samtoolsModule + " stats " + bamFilePath + " >> " + statsOutFile + "\n")
    shellFile.close()
    os.chmod(shScriptPath, 0o744)

def makeJobListFile():
    jobListFile = open(jobFilePath, "w")
    for i in shScriptPathDict:
        jobListFile.write(shScriptPathDict[i] + "\n")
    jobListFile.close()

# Generate a SLURM file for the scripts mentioned in the JobListFile
def makepbsJobFile():
    pbsFile = open(pbsDirName + "/reads_extract.pbs", "w")
    pbsFile.write("#!/bin/bash\n\n")
    pbsFile.write("#PBS -V\n")
    pbsFile.write("#PBS -N reads_extract\n")
    pbsFile.write("#PBS -q default\n")    
    pbsFile.write("#PBS -l walltime=200:00:00\n")
    pbsFile.write("#PBS -o " + logPath + "/reads_extract_pbs.log\n")
    pbsFile.write("#PBS -l nodes=1:ppn=20\n")
    pbsFile.write("#PBS -M jcvalverdehernandez@gmail.com\n")
    pbsFile.write("#PBS -m abe\n\n")
    pbsFile.write("export OMP_NUM_THREADS=10\n")

    pbsFile.write("parallel=\"parallel --delay .2 -j 10 --joblog " + logPath + "/reads_extract_parallel.log \"\n")
    pbsFile.write("$parallel < " + jobFilePath + "\n")
    pbsFile.close()
    pbsFile.close()

### IMPLEMENTATION


#Create shell scripts
jobFilePath = joblistsDirName + "/JobList.txt"
for i in cramDirFiles:
    createScripts(i)

# Create job list
makeJobListFile()

# Create slurm file
makepbsJobFile()

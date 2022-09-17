#!/usr/bin/python


### IMPORT

import os
import os.path
import sys

### GLOBAL VARIABLES

# Location of the path that will contain the scripts for this phase
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/
scriptsDir = sys.argv[1]
# Location of folder containing the BAM files to be genotyped and unified
# i.e. /home/jcvh/cr_dislipidemia/data/03_recalibrateBAM/
bamDir = sys.argv[2]
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
# and index (fai) need to be stored
# i.e. /home/jcvh/cr_dislipidemia/study-specs/genomicCoordinates/wg.list
coordinatesListPath = sys.argv[6]
# Name of the GATK module
# i.e. /opt/gatk-4.2.0.0/gatk
gatkPath = sys.argv[7]
# Name of java module
# i.e. jdk/13.0.2
##javaModule = sys.argv[8]

# Identifier of the samples in the bamDir folder
bamDirFilesTemp = os.listdir(bamDir)
bamDirFiles = []
for i in bamDirFilesTemp:
    file_name = i.split('.')
    if file_name[-1] == "bam":
        bamDirFiles.append('.'.join(file_name))
    else:
        pass
# Location to save generated pbs script
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/shScripts/
shScriptDirName = scriptsDir + "shScripts/"
# Location to save generated scripts list
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/jobLists
joblistsDirName = scriptsDir + "jobLists/"
# Location to save generated pbs script
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/pbs/
pbsDirName = scriptsDir + "pbs/"
# Location of the GVCF file map
# i.e. /home/jcvh/cr_dislipidemia/logs/04_createGenotype/gvcfMap.tsv
gvcfMapPath = logPath + "gvcfMap.tsv"
#Path to sh job list
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/jobLists/jobList.txt
jobFilePath = joblistsDirName + "jobList.txt"

##MAKE DIRECTORIES
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/shScripts/
if (not os.path.isdir(shScriptDirName)):
    os.mkdir(shScriptDirName)
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/jobLists/
if (not os.path.isdir(joblistsDirName)):
    os.mkdir(joblistsDirName)
# i.e. /home/jcvh/cr_dislipidemia/scripts/04_createGenotype/pbs/
if (not os.path.isdir(pbsDirName)):
    os.mkdir(pbsDirName)


### FUNCTION DEFINITIONS

gvcfPointersDict = {}
def createScripts(bamID):
    # Create the shell script file
    bamFilePath = bamDir + bamID
    bamID = bamID.split('.')[0]
    shScriptPath = shScriptDirName + bamID + ".sh"
    individualGVCFFilePath = dataSavePath + bamID + ".gvcf.gz"
    gvcfPointersDict[bamID] = [shScriptPath, individualGVCFFilePath]
    haploCallLogPath = logPath + bamID + "_haploCall.log"
    sha256sumLogPath = logPath + "sha256sum.log"
    md5sumLogPath = logPath + "md5sum.log"

    shellFile = open(shScriptPath, "w")
    # Write to the file
    shellFile.write("#!/bin/bash\n\n")
    shellFile.write("start=`date +%s`\n")
    shellFile.write("start_moment=$(date)\n")
    shellFile.write(gatkPath + ' --java-options "-Xmx4g" HaplotypeCaller -ERC GVCF -R ' + refPath + " -I " + bamFilePath + " -O " + individualGVCFFilePath + " -G StandardAnnotation -G StandardHCAnnotation &> " + haploCallLogPath + "\n")
    shellFile.write("sha256sum " + individualGVCFFilePath + " >> " + sha256sumLogPath + "\n")
    shellFile.write("md5sum " + individualGVCFFilePath + " >> " + md5sumLogPath + "\n")
    shellFile.write("end=`date +%s`\n")
    shellFile.write("end_moment=$(date)\n")
    shellFile.write("echo Job for " + bamID + " started on $start_moment and ended at $end_moment - Total time: $((end-start))\n")
    shellFile.close()
    os.chmod(shScriptPath, 0o744)

def makeJobListFile():
    jobListFile = open(jobFilePath, "w")
    for i in gvcfPointersDict:
        jobListFile.write(gvcfPointersDict[i][0] + "\n")
    jobListFile.close()
    
# Generate a pbs file for the scripts mentioned in the JobListFile
def makepbsJobFile():
    pbsFile = open(pbsDirName + "/variants_individual.pbs", "w")
    pbsFile.write("#!/bin/bash\n\n")
    pbsFile.write("#PBS -V\n")
    pbsFile.write("#PBS -N variants_individual\n")
    pbsFile.write("#PBS -q default\n")    
    pbsFile.write("#PBS -l walltime=200:00:00\n")
    pbsFile.write("#PBS -o " + logPath + "/variants_individual_pbs.log\n")
    #pbsFile.write("#PBS --q dribe\n") all partitions are the same
    pbsFile.write("#PBS -l nodes=1:ppn=" + "20" + "\n")
#    pbsFile.write("#PBS --exclusive\n")
    pbsFile.write("#PBS -M jcvalverdehernandez@gmail.com\n")
    pbsFile.write("#PBS -m abe\n\n")
    pbsFile.write("export OMP_NUM_THREADS=6\n")

    pbsFile.write("echo Start time: $(date)\n")
    pbsFile.write("parallel=\"parallel --delay .2 -j 6 --joblog " + logPath + "/variants_individual_parallel.log --resume\"\n")
    pbsFile.write("$parallel < " + jobFilePath + "\n")
    pbsFile.write("echo End time: $(date)\n")
    pbsFile.close()


### IMPLEMENTATION

for i in bamDirFiles:
    createScripts(i)

# Create job list
makeJobListFile()

# Create pbs file
makepbsJobFile()

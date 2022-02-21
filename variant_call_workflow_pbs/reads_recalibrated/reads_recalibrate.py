#!/usr/bin/python


### IMPORT

import os
import os.path
import sys

### GLOBAL VARIABLES

# Location of the path that will contain the scripts for this phase
# i.e. /home/jcvh/cr_dislipidemia/scripts/03_recalibrateBAM/
scriptsDir = sys.argv[1]
# Location of folder containing the BAM files to recalibrate - focused for ILLUMINA
# i.e. /home/jcvh/cr_dislipidemia/data/02_createBAM/
bamDir = sys.argv[2]
# Path to store log and intermediate files
# i.e. /home/jcvh/cr_dislipidemia/logs/03_recalibrateBAM/
logPath = sys.argv[3]
# Storage location of recalibrated BAM files
# i.e. /home/jcvh/cr_dislipidemia/data/03_recalibrateBAM/
dataSavePath = sys.argv[4]
# Path to the reference genome in format. In the same directory its dictionary
# and index (fai) need to be stored
# i.e. /home/jcvh/cr_dislipidemia/study-specs/referenceGenome/Homo_sapiens_assembly38.fasta 
refPath = sys.argv[5]
# Path to the VCF file containing the known variant sites (such as dbSNP or 1k genomes database).
# Its respective index file needs to be stored in the same folder (generated with gatk IndexFeatureFile)
# i.e. /home/jcvh/cr_dislipidemia/study-specs/variantsReference/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
variantsPath = sys.argv[6]
# Location of the GATK jar
# i.e.  /opt/gatk-4.2.0.0/gatk
gatkPath = sys.argv[7]


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
# i.e. /work/jvalverde/cr_dislipidemia/scripts/03_recalibrateBAM/shScripts/
shScriptDirName = scriptsDir + "shScripts/"
# Location to save generated scripts list
# i.e. /work/jvalverde/cr_dislipidemia/scripts/03_recalibrateBAM/jobLists
joblistsDirName = scriptsDir + "jobLists/"
# Location to save generated pbs script
# i.e. /work/jvalverde/cr_dislipidemia/scripts/03_recalibrateBAM/pbs/
pbsDirName = scriptsDir + "pbs/"


##MAKE DIRECTORIES
# i.e. /work/jvalverde/cr_dislipidemia/scripts/03_recalibrateBAM/shScripts/
if (not os.path.isdir(shScriptDirName)):
    os.mkdir(shScriptDirName)
# i.e. /work/jvalverde/cr_dislipidemia/scripts/03_recalibrateBAM/jobLists/
if (not os.path.isdir(joblistsDirName)):
    os.mkdir(joblistsDirName)
# i.e. /work/jvalverde/cr_dislipidemia/scripts/03_recalibrateBAM/pbs/
if (not os.path.isdir(pbsDirName)):
    os.mkdir(pbsDirName)

### FUNCTION DEFINITIONS

shScriptPathList = []
def createScripts(bamID):
    # Create the shell script file
    bamFilePath = bamDir + bamID
    bamID = bamID.split('.')[0]
    shScriptPath = shScriptDirName + bamID + ".sh"
    shScriptPathList.append(shScriptPath)
    bamFilePathNoDuplicates = bamDir + bamID + "_nodup.bam"
    MarkDuplicatesuplicatesMetrics = logPath + bamID + "_DuplicateMetrics.txt"
    MarkDuplicatesLogFile = logPath + bamID + "_MarkDuplicates.log"
    baseRecalTablePath_1 = logPath + bamID + "_baseRecal1.table"
    baseRecalLogPath_1 = logPath + bamID + "_baseRecal1.log"
    bamRecalibratedPath = dataSavePath + bamID + ".recal.bam"
    applyBQSRLogPath = logPath + bamID + "_applyBQSR.log"
    baseRecalTablePath_2 = logPath + bamID + "_baseRecal2.table"
    baseRecalLogPath_2 = logPath + bamID + "_baseRecal2.log"
    covariantsCsvPath =logPath + bamID + "_analyzeCovariants.csv"
    covariantsLogPath =logPath + bamID + "_analyzeCovariants.log"
    sha256sumLogPath = logPath + "sha256sum.log"
    md5sumLogPath = logPath + "md5sum.log"
    shellFile = open(shScriptPath, "w")
    # Write to the file
    shellFile.write("#!/bin/bash\n\n")
    shellFile.write("start=`date +%s`\n")    
    shellFile.write(gatkPath + " MarkDuplicates -I " + bamFilePath + " -O " + bamFilePathNoDuplicates + " -M " + MarkDuplicatesuplicatesMetrics + " --REMOVE_DUPLICATES true &> " + MarkDuplicatesLogFile + "\n")
    shellFile.write(gatkPath + " BaseRecalibrator -I " + bamFilePathNoDuplicates + " -R " + refPath + " --known-sites " + variantsPath + " -O " + baseRecalTablePath_1 + " &> " + baseRecalLogPath_1 + "\n")
    shellFile.write(gatkPath + " ApplyBQSR -I " + bamFilePathNoDuplicates + " -R " + refPath + " --bqsr-recal-file " + baseRecalTablePath_1 + " -O " + bamRecalibratedPath + " &> " + applyBQSRLogPath + "\n")
    shellFile.write(gatkPath + " BaseRecalibrator -I " + bamRecalibratedPath + " -R " + refPath + " --known-sites " + variantsPath + " -O " + baseRecalTablePath_2 + " &> " + baseRecalLogPath_2 + "\n")
    shellFile.write(gatkPath + " AnalyzeCovariates -before " + baseRecalTablePath_1 + " -after " + baseRecalTablePath_2 + " -csv " + covariantsCsvPath + " &> " + covariantsLogPath + "\n")
    shellFile.write("sha256sum " + bamRecalibratedPath + " >> " + sha256sumLogPath + "\n")
    shellFile.write("md5sum " + bamRecalibratedPath + " >> " + md5sumLogPath + "\n")
    shellFile.write("end=`date +%s`\n")
    shellFile.write("echo " + bamID + "\n")    
    shellFile.write("echo $((end-start))\n")   
    shellFile.close()
    os.chmod(shScriptPath, 0o744)

def makeJobListFile():
    jobListFile = open(jobFilePath, "w")
    for i in shScriptPathList:
        jobListFile.write(  i + "\n")
    jobListFile.close()

# Generate a pbs file for the scripts mentioned in the JobListFile
def makepbsJobFile():
    pbsFile = open(pbsDirName + "/reads_recalibrated.pbs", "w")
    pbsFile.write("#!/bin/bash\n\n")
    pbsFile.write("#PBS -V\n")
    pbsFile.write("#PBS -N reads_recalibrated\n")
    pbsFile.write("#PBS -q default\n")    
    pbsFile.write("#PBS -l walltime=200:00:00\n")
    pbsFile.write("#PBS -o " + logPath + "/reads_recalibrated_pbs.log\n")
    pbsFile.write("#PBS -l nodes=1:ppn=20\n")
#    pbsFile.write("#PBS --exclusive\n")
    pbsFile.write("#PBS -M jcvalverdehernandez@gmail.com\n")
    pbsFile.write("#PBS -m abe\n\n")
    pbsFile.write("export OMP_NUM_THREADS=4\n")

    pbsFile.write("parallel=\"parallel --delay .2 -j 4 --joblog " + logPath + "/reads_recalibrated_parallel.log --resume\"\n")
    pbsFile.write("$parallel < " + jobFilePath + "\n")
    pbsFile.close()

### IMPLEMENTATION

#Create shell scripts
jobFilePath = joblistsDirName + "jobList.txt"
for i in bamDirFiles:
    createScripts(i)

# Create job list
makeJobListFile()

# Create pbs file
makepbsJobFile()

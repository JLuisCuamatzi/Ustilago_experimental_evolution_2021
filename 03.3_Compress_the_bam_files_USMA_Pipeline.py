###
###
###############################################################################################################################################
### File:           03.3_Compress_the_bam_files_USMA_Pipeline.py
### Written by:     Jorge Luis Cuamatzi Flores
### Date:           2020_May_29
### Update:         2021_September_19
###
### Project:        This script writes a sge file to map and do the variant calling
### Input:          A csv file with information and simple paths
### Output:         SGE to do mapping and variant calling for several samples
###
# How execute this script: python3 /mnt/Timina/lmorales/Public/Ustilago/C1/bin/scripts/03.3_Compress_the_bam_files_USMA_Pipeline.py -d /mnt/Timina/lmorales/Public/Ustilago/C1/ -f /mnt/Timina/lmorales/Public/Ustilago/C1/ID.csv -t Mapping -r /mnt/Timina/lmorales/Public/Ustilago/reference/USMA_521_v2.csv -M 8
###############################################################################################################################################
## Libraries
import argparse
import os
from datetime import date, datetime
import sys
import subprocess
##
fecha = date.today().strftime('%Y%m%d')
tiempo = datetime.now().strftime("%H.%M")
##
ag = argparse.ArgumentParser()
ag = argparse.ArgumentParser(
    description = "Python scripts that takes information from csv file and write lines",
    usage = "python3 Practice_script.py -f ~/USMA_ExpEvo_Samples.csv")
ag.add_argument("-f", "--file", default = "", help = "csv with information of ID") # to read a csv file with sample imformation
ag.add_argument("-n", "--email", default = "jorge_l.1@hotmail.com", help = "If you want to receive a notification when the process is done") 
ag.add_argument("-r", "--ref", default = "", help = "csv file with the data of the references to map to")
ag.add_argument("-d", "--directory", default = "/home/jcuamatzi", help = "path to the project directory")
ag.add_argument("-t", "--task", default = "", help = "Task")
ag.add_argument("-M", "--memory", default = "2", help = "RAM memory")

##
#
args = vars(ag.parse_args())
arg_file = args["file"]
arg_reference = args["ref"]
email = str(args["email"])
directory = str(args["directory"])
task = str(args["task"])
memory = str(args["memory"])

##### Functions #####
# Function to open a csv file
# With this function, we are only open the file
def opencsv (csv_file):
    arr = []
    csv_file = open(csv_file,"r",encoding = "utf-8")
    for line in csv_file:
        arr.append(line.strip().split(","))
    csv_file.close
    return arr
# Function to look at info from the header
# This is a generic function in order to look at the header (first row) of any csv file
def extindex (row,header):
    col_ID = 0
    for col in row:
        if header in col:
            break
        else:
            col_ID += 1        
    return col_ID
# Function to extract info from the column at a specific header
# This function employ the extindex function
def extcol (array, header):
    index = extindex(array[0], header)
    arr = []
    for i in range(1, len(array)):
        for j in range(0, len(array[0])):
            if j == index:
                arr.append(array[i][j])
            else:
                pass
    return arr
# Function to write the header of sge file
#
#def header(smpls_ID, sge):
def header(smpls_ID,sge):
    path_error = wd_project + "log/" + task + "/error/"
    path_out = wd_project + "log/" + task + "/out/"
    if not os.path.exists(path_error):
        os.makedirs(path_error)
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    print('''#!/bin/bash
## Use current working directory
#$ -cwd
. /etc/profile.d/modules.sh
##Error file
#$ -e''', path_error + smpls_ID + "." + fecha + "." + tiempo + ".error",'''
##Out file
#$ -o''', path_out + smpls_ID + "." + fecha + "." + tiempo + ".out",'''
#$ -S /bin/bash
## Job's name
#$ -N''', task + "_" + ID,'''
#$ -l vf='''+ memory +'''G
#$ -pe openmp 10
#$ -m e
source /etc/bashrc
## notification
#$ -M ''' + email + '''
##
## Modules''', file = sge)
    print ("module load htslib/1.2.1 gcc/5.1.0 samtools/1.9", file = sge)
    print ("##", file = sge)
    print ("##", file = sge)

def mapping (smpls_ID):
    save_sge = wd_project + "bin/SGE/" + task + "/Compress_bam_files"
    if not os.path.exists(save_sge):    
        os.makedirs(save_sge)
    sge_name = save_sge + "/" + fecha + "_" + smpls_ID + "_" + sample_name + "_compress_bam_files.sge"
    sge = open(sge_name, "w")
    header(smpls_ID,sge)
    bam_path = wd_project + "data/bam/"
    ref_file = wd_ref + "/" + ref_name + ".fa"
    # Names
    bam_name = smpls_ID + "_" + sample_name + "_BWA.bam"
    ## Coverage requirements
    # Coverage requirements: paths
    #
    print ('''start=$(date +%s.%N)''', file = sge)
    ## Coverage
    print ("## START Compress the bam files", file = sge)
    print ("samtools view -T " + ref_file + " -C -o " + bam_path + smpls_ID + ".cram" + " " + bam_path + bam_name, file = sge)
    print ("rm " + bam_path + bam_name, file = sge)
    print ("## END Compress the bam files", file = sge)
    print ('''duration=$(echo "$(date +%s.%N) - $start" | bc)''', file = sge)
    print ('''execution_time=`printf "%.2f seconds" $duration`''', file = sge)
    print ('''echo "Script Execution Time: $execution_time"''', file = sge)
    sge.close()
    return sge_name

##### MAIN #####
## OPEN FILE ##
matrix_csv = opencsv(arg_file)
reference_csv = opencsv(arg_reference)
wd_project = directory
## IDs & NAMES ##
smpls_ID = extcol(matrix_csv, "#ID")
wd_ref = extcol(reference_csv, "RefPath")[0]
pyout_name = wd_project + "/log/Py_scripts/" + task + "/" + fecha + "_" + task + "_py.out"
if not os.path.exists(wd_project + "/log/Py_scripts/" + task + "/"):
            os.makedirs(wd_project + "/log/Py_scripts/" + task + "/")
pyout_name = os.path.normpath(pyout_name)
pyoutput = open(pyout_name, "a+")
for i in range(0, len(smpls_ID)):
    ID = smpls_ID[i]
    sample_name = extcol(matrix_csv, "Name")[i]
    ref_name = extcol(reference_csv, "RefName")[0]
    sge = mapping(ID)
    print(tiempo, file = pyoutput)
    #subprocess.run(["qsub",sge], stdout=pyoutput)
    print(tiempo, file = pyoutput)
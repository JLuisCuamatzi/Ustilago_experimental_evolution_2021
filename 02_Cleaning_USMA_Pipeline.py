###
###
###############################################################################################################################################
### File:           02_Cleaning_USMA_Pipeline.py
### Written by:     Jorge Luis Cuamatzi Flores
### Date:           2021_August_13
###
### Project:        This script writes a sge file to map and do the variant calling
### Input:          A csv file with information and simple paths
### Output:         SGE to do mapping and variant calling for several samples
###
# How execute this script: python3 /mnt/Timina/lmorales/Public/Ustilago/C1/bin/scripts/02_Cleaning_USMA_Pipeline.py -d /mnt/Timina/lmorales/Public/Ustilago/C1/ -f /mnt/Timina/lmorales/Public/Ustilago/C1/ID.csv -t Cleaning -w 10 -M 6 -n jcuamatzi@liigh.unam.mx
# -w threads for fastp
# -M RAM memory
# for the -t cannot start with a digit
###############################################################################################################################################
## Libraries
import argparse
import os
from datetime import date, datetime
import sys
import subprocess
##
fecha = date.today().strftime('%Y%m%d')
tiempo = datetime.now().strftime("%H:%M:%S")
##
ag = argparse.ArgumentParser()
ag = argparse.ArgumentParser(
    description = "Python scripts that takes information from csv file and write lines",
    usage = "python3 Practice_script.py -f ~/USMA_ExpEvo_Samples.csv")
ag.add_argument("-f", "--file", default = "", help = "csv with information of ID") # to read a csv file with sample imformation
ag.add_argument("-n", "--email", default = "jcuamatzi@liigh.unam.mx", help = "If you want to receive a notification when the process is done") 
ag.add_argument("-d", "--directory", default = "/home/jcuamatzi/", help = "path to the project directory")
ag.add_argument("-t", "--task", default = "", help = "Task")
ag.add_argument("-w", "--thread", default = "2", help = "thread number, default is 2")
ag.add_argument("-M", "--memory", default = "2", help = "RAM memory")
##
#
args = vars(ag.parse_args())
arg_file = args["file"]
email = str(args["email"])
directory = str(args["directory"])
task = str(args["task"])
thread = str(args["thread"])
memory = str(args["memory"]) #RAM Memory
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
#$ -e''', path_error + ID + ".error",'''
##Out file
#$ -o''', path_out + smpls_ID + ".out",'''
#$ -S /bin/bash
## Job's name
#$ -N''', task + "_" + ID,'''
#$ -l vf='''+ memory +'''G
#$ -pe openmp 2
#$ -m e
source /etc/bashrc
## notification
#$ -M ''' + email + '''
##
## Modules''', file = sge)
    print ("module load fastp/0.20.0", file = sge)
    print ("##", file = sge)
    print ("##", file = sge)

def cleaning(smpls_ID):
    save_sge = wd_project + "bin/SGE/" + task + "/"
    if not os.path.exists(save_sge):    
        os.makedirs(save_sge)
    sge_name = save_sge + fecha + "_" + smpls_ID + "_" + sample_name + ".sge"
    sge = open(sge_name, "w")
    header(smpls_ID,sge)
    ## Cleaning requirements
    ## Paths
    fastq_raw_path = wd_project + "data/fastq/raw/"
    fastq_clean_path = wd_project + "data/fastq/clean/"
    if not os.path.exists(fastq_clean_path):
        os.makedirs(fastq_clean_path)
    fastq_unpaired_path = wd_project + "data/fastq/clean/unpaired/"
    if not os.path.exists(fastq_unpaired_path):
        os.makedirs(fastq_unpaired_path)
    fastp_json_path = wd_project + "tmp/" + task + "/"
    if not os.path.exists(fastp_json_path):
        os.makedirs(fastp_json_path)
    fastp_htlml_path = wd_project + "tmp/" + task + "/"
    if not os.path.exists(fastp_htlml_path):
        os.makedirs(fastp_htlml_path)
    crr_path = wd_project + "analysis/reads_num/raw/" #count raw reads
    if not os.path.exists(crr_path):
        os.makedirs(crr_path)
    ccr_path = wd_project + "analysis/reads_num/clean/" #count clean reads
    if not os.path.exists(ccr_path):
        os.makedirs(ccr_path)
    # Names
    fastq_name_raw_R1 = smpls_ID + "_1.fq.gz" # for compress fastq
    #fastq_name_raw_R1 = smpls_ID + "_" + sample_name + "_R1.fastq" # for uncompress fastq
    fastq_name_raw_R2 = smpls_ID + "_2.fq.gz" # for compress fastq
    #fastq_name_raw_R2 = smpls_ID + "_" + sample_name + "_R2.fastq" # for uncompress fastq
    fastq_name_clean_R1 = smpls_ID + "_" + sample_name + "_R1_clean.fq.gz"
    fastq_name_clean_R2 = smpls_ID + "_" + sample_name + "_R2_clean.fq.gz"
    fastq_name_unpaired = smpls_ID + "_" + sample_name + "_unpaired_clean.fastq.gz"
    json_name = smpls_ID + "_fastp.json"
    html_name = smpls_ID + "_fastp.html"
    #
    print ('''start=$(date +%s.%N)''', file = sge)
    print ("# Calculating the number of reads in the raw fastq file",file = sge)
    print ("echo $(zcat " + fastq_raw_path + fastq_name_raw_R1 + " | wc -l )/4 | bc > " + crr_path + smpls_ID + "_raw_reads.txt",file = sge)
    print ("echo $(zcat " + fastq_raw_path + fastq_name_raw_R2 + " | wc -l )/4 | bc > " + crr_path + smpls_ID + "_raw_reads.txt",file = sge)
    print ("## CLEANING ", file = sge)
    print ('''echo "Start to clean fastq files"''', file = sge)
    print ("#", file = sge)
    print ("fastp -i " + fastq_raw_path + fastq_name_raw_R1 + " -o " + fastq_clean_path + fastq_name_clean_R1 + " -I " + fastq_raw_path + fastq_name_raw_R2 + " -O " + fastq_clean_path + fastq_name_clean_R2 + " --unpaired1 " + fastq_unpaired_path + fastq_name_unpaired + " -w " + thread + " -y -x -z 9 -j " + fastp_json_path + json_name + " -h " + fastp_htlml_path + html_name, file = sge)
    print ("#", file = sge)
    print ("#Calculating the number of reads in the clean fastq file", file = sge)
    print ("echo $(zcat " + fastq_clean_path + fastq_name_clean_R1 + " | wc -l )/4 | bc > " + ccr_path + smpls_ID + "_clean_reads.txt",file = sge)
    print ("echo $(zcat " + fastq_clean_path + fastq_name_clean_R2 + " | wc -l )/4 | bc > " + ccr_path + smpls_ID + "_clean_reads.txt",file = sge)
    print ('''duration=$(echo "$(date +%s.%N) - $start" | bc)''', file = sge)
    print ('''execution_time=`printf "%.2f seconds" $duration`''', file = sge)
    print ('''echo "Script Execution Time: $execution_time"''', file = sge)
    print ('''echo "The script ends"''', file=sge)
    sge.close()
    return sge_name  

##### MAIN #####
## OPEN FILE ##
matrix_csv = opencsv(arg_file)
wd_project = directory
## IDs & NAMES ##
smpls_ID = extcol(matrix_csv, "#ID")
pyout_name = wd_project + "/log/Py_scripts/" + task + "/" + fecha + "_" + task + "_py.out"
if not os.path.exists(wd_project + "/log/Py_scripts/" + task + "/"):
            os.makedirs(wd_project + "/log/Py_scripts/" + task + "/")
pyout_name = os.path.normpath(pyout_name)
pyoutput = open(pyout_name, "a+")
for i in range(0, len(smpls_ID)):
    ID = smpls_ID[i]
    sample_name = extcol(matrix_csv, "Name")[i]
    sge = cleaning(ID)
    print(tiempo, file = pyoutput)
    #subprocess.run(["qsub",sge], stdout=pyoutput)
    print(tiempo, file = pyoutput)
###
###
###############################################################################################################################################
### File:           01_Downloading_fq_files_USMA_Pipeline.py
### Written by:     Jorge Luis Cuamatzi Flores
### Date:           2021_September_16
###
### Project:        This script writes a sge file to download the fastq files from BGI cloud
### Input:          A csv file with information and simple paths
### Output:         SGE to do mapping and variant calling for several samples
###
# How execute this script: python3 /mnt/Timina/lmorales/Public/Ustilago/C1/bin/scripts/01_Downloading_fq_files_USMA_Pipeline.py -d /mnt/Timina/lmorales/Public/Ustilago/C1/ -f /mnt/Timina/lmorales/Public/Ustilago/C1/ID.csv -t DWN -m 6 -n jcuamatzi@liigh.unam.mx -s s3://funmnljr/F20FTSUSAT1199_FUNmnljR/Filter_SOAPnuke/Clean
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
ag.add_argument("-m", "--memory", default = "2", help = "RAM memory")
ag.add_argument("-s", "--aws", default = "", help = "path for the AWS cloud")
##
#
args = vars(ag.parse_args())
arg_file = args["file"]
email = str(args["email"])
directory = str(args["directory"])
task = str(args["task"])
memory = str(args["memory"]) #RAM Memory
aws_path = str(args["aws"]) # aws path
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
    print ("module load awscli/2.0.9", file = sge)
    print ("##", file = sge)
    print ("##", file = sge)

def downloading(smpls_ID):
    save_sge = wd_project + "bin/SGE/" + task + "/"
    if not os.path.exists(save_sge):    
        os.makedirs(save_sge)
    sge_name = save_sge + fecha + "_" + smpls_ID + "_" + sample_name + ".sge"
    sge = open(sge_name, "w")
    header(smpls_ID,sge)
    ## Cleaning requirements
    ## Paths
    fastq_raw_path = wd_project + "data/fastq/raw/" + ID + "/"
    if not os.path.exists(fastq_raw_path):
        os.makedirs(fastq_raw_path)
    fq_cloud_path = aws_path + "/" + ID
    # Names
    fastq_name_raw_R1 = smpls_ID + "_1.fq.gz" # for compress fastq
    #fastq_name_raw_R1 = smpls_ID + "_" + sample_name + "_R1.fastq" # for uncompress fastq
    fastq_name_raw_R2 = smpls_ID + "_2.fq.gz" # for compress fastq
    #fastq_name_raw_R2 = smpls_ID + "_" + sample_name + "_R2.fastq" # for uncompress fastq
    #
    print ("## DOWNLOADING ", file = sge)
    print ('''start=$(date +%s.%N)''', file = sge)
    print ('''echo "Start to download the fastq files"''', file = sge)
    print ("#", file = sge)
    print ("aws s3 sync " + fq_cloud_path + " " + fastq_raw_path, file = sge)
    print ("#", file = sge)
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
    sge = downloading(ID)
    print(tiempo, file = pyoutput)
    subprocess.run(["qsub",sge], stdout=pyoutput)
    print(tiempo, file = pyoutput)
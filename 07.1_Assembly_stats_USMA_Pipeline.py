###
###
###############################################################################################################################################
### File:           07.1_Assembly_stats_USMA_Pipeline.py
### Written by:     Jorge Luis Cuamatzi Flores
### Date:           2021_October_14
### Update:         
###
### Project:        This script writes a sge file to do a Copy Number Variation with CNVator
### Input:          A csv file with information and simple paths
### Output:         SGE to do mapping and variant calling for several samples
###
# How execute this script:
#python3 /mnt/Timina/lmorales/Public/Ustilago/C1/bin/scripts/07.1_Assembly_stats_USMA_Pipeline.py -d /mnt/Timina/lmorales/Public/Ustilago/C1/ -f /mnt/Timina/lmorales/Public/Ustilago/C1/ID.csv -t Assembly -M 8 -f /mnt/Timina/lmorales/Public/Ustilago/C1/ID.csv -p /mnt/Timina/lmorales/Public/Ustilago/C1/data/Assembly/
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
ag.add_argument("-t", "--task", default = "", help = "Task of this script")
ag.add_argument("-d", "--directory", default = "/home/jcuamatzi", help = "path to the project directory")
ag.add_argument("-M", "--memory", default = "2", help = "RAM memory")
ag.add_argument("-p", "--pathAssembly", default = "", help = "path for the fastq files")
##
#
args = vars(ag.parse_args())
arg_file = args["file"]
email = str(args["email"])
task = str(args["task"])
directory = str(args["directory"])
memory = str(args["memory"]) #RAM Memory
path_assembly = str(args["pathAssembly"])
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
#def header(GeneByChr, sge):
def header(smpls_ID, sge):
    path_error = wd_project + "log/" + task + "/error/Stats/"
    path_out = wd_project + "log/" + task + "/out/Stats/"
    if not os.path.exists(path_error):
        os.makedirs(path_error)
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    print('''#!/bin/bash
## Use current working directory
#$ -cwd
. /etc/profile.d/modules.sh
##Error file
#$ -e''', path_error + ID + "_" + task + "_Stats.error",'''
##Out file
#$ -o''', path_out + ID + "_" + task + "_Stats.out",'''
#$ -S /bin/bash
## Job's name
#$ -N''', task + "_" + smpls_ID,'''
#$ -l vf='''+ memory +'''G
#$ -pe openmp 10
#$ -m e
source /etc/bashrc
## notification
#$ -M ''' + email + '''
##''', file = sge)
    print ("module load assembly-stats/1.0.1", file = sge)
    print ("##", file = sge)

def mapping (smpls_ID):
    save_sge = wd_project + "bin/SGE/" + task + "/Stats"
    if not os.path.exists(save_sge):    
        os.makedirs(save_sge)
    sge_name = save_sge + "/" + fecha + "_" + smpls_ID + "_DeNovo_Assembly_Stats.sge"
    sge = open(sge_name, "w")
    header(smpls_ID, sge)
    # paths
    output_path = wd_project + "analysis/DeNovoAssembly/"
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    print ('''start=$(date +%s.%N)''', file = sge)
    print ("## START: Genome Assembly Stats", file = sge)
    print ("mv " + path_assembly + smpls_ID + "/scaffolds.fasta " + path_assembly + smpls_ID + "/" + smpls_ID + "_scaffolds.fasta", file = sge)
    print ("rm " + path_assembly + smpls_ID + "/scaffolds.fasta", file = sge)
    print ("assembly-stats " + path_assembly + smpls_ID + "/" + smpls_ID + "_scaffolds.fasta > " + output_path + smpls_ID + "_stats_DeNovoAssembly.txt", file = sge)
    print ("## END: Genome Assembly Stats", file = sge)
    print ('''duration=$(echo "$(date +%s.%N) - $start" | bc)''', file = sge)
    print ('''execution_time=`printf "%.2f seconds" $duration`''', file = sge)
    print ('''echo "Script Execution Time: $execution_time"''', file = sge)
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
    sge = mapping(ID)
    print(tiempo, file = pyoutput)
    subprocess.run(["qsub",sge], stdout=pyoutput)
    print(tiempo, file = pyoutput)
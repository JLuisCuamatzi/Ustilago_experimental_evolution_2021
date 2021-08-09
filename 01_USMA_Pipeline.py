###
###
###############################################################################################################################################
### File:           01_USMA_Pipeline.py
### Written by:     Jorge Luis Cuamatzi Flores
### Date:           2020_May_29
### Update:         2021_July_27
###
### Project:        This script writes a sge file to map and do the variant calling
### Input:          A csv file with information and simple paths
### Output:         SGE to do mapping and variant calling for several samples
###
# How execute this script: python3 /mnt/Cromosoma/lmorales/Public/Ustilago/B1/bin/scripts/python_scripts/01_JLCF_Mapping.py -f /mnt/Cromosoma/lmorales/Public/Ustilago/B1/information/USMA_ExpEvo_Samples.csv -r /mnt/Cromosoma/lmorales/Public/Ustilago/B1/information/ref_USMA_521_v2.csv -n jcuamatzi@dna.lavis.unam.mx -m BWA
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
ag.add_argument("-n", "--email", default = "jorge_l.1@hotmail.com", help = "If you want to receive a notification when the process is done") 
ag.add_argument("-r", "--ref", default = "", help = "csv file with the data of the references to map to")
ag.add_argument("-m", "--mapping", default = "BWA", help = "mapping tool. BWA or Bowtie2 (B2)")
# agregar un ardgumento para el working directory
ag.add_argument("-d", "--directory", default = "/home/jcuamatzi", help = "path to the project directory")
ag.add_argument("-t", "--task", default = "", help = "Task")
##
#
args = vars(ag.parse_args())
arg_file = args["file"]
arg_reference = args["ref"]
email = str(args["email"])
map_tool = str(args["mapping"])
directory = str(args["directory"])
task = str(args["task"])
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
#$ -e''', path_error + "/" + ID + ".error",'''
##Out file
#$ -o''', path_out  + "/" + smpls_ID + ".out",'''
#$ -S /bin/bash
## Job's name
#$ -N''', task + "_" + ID,'''
#$ -l vf=8G
#$ -pe openmp 10
#$ -m e
source /etc/bashrc
## notification
#$ -M ''' + email + '''
##
## Modules''', file = sge)
    if map_tool == "BWA" or map_tool == "bwa":
        print ("module load bwa/0.7.4 htslib/1.2.1 gcc/5.1.0 samtools/1.9 picard/2.6.0 r/3.6.1 bedops/2.4.20 bedtools/2.24 gatk/4.1.1.0", file = sge)
    elif map_tool == "B2" or map_tool == "Bowtie2":
        print ("module load htslib/1.2.1 gcc/5.1.0 samtools/1.9 python37/3.7.0 fastqc/0.11.3 picard/2.6.0 bamtools/2.5.1 bowtie2/2.2.9", file = sge)
    else:
        sys.exit("No mapping tool selected, bye-bye")
    print ("##", file = sge)
    print ("##", file = sge)

def pipeline(smpls_ID):
    save_sge = wd_project + "/bin/SGE/" + task + "/"
    if not os.path.exists(save_sge):    
        os.makedirs(save_sge)
    sge_name = save_sge + "/" + fecha + "_" + smpls_ID + ".sge"
    sge = open(sge_name, "w")
    header(smpls_ID,sge)
    #reference_genome = wd_ref + "/" + ref_name + ".fa"
    ## Mapping requirements
    ## Paths
    #bam_path = "/mnt/Cromosoma/lmorales/Public/Ustilago/B1/data/bam"
    #fastq_path = "/mnt/Cromosoma/lmorales/Public/Ustilago/A1/data/fastq/clean/good"
    #bam_stats_path = bam_path + "/stats"
    #if not os.path.exists(bam_stats_path):
    #    os.makedirs(bam_stats_path)
    #stats_path = "/mnt/Cromosoma/lmorales/Public/Ustilago/B1/analysis/stats"
    #mapp_stats_path = stats_path + "/mapp"
    #if not os.path.exists(mapp_stats_path):
    #    os.makedirs(mapp_stats_path)
    #stats_GC_path = mapp_stats_path + "/GCBias"
    #if not os.path.exists(stats_GC_path):
    #    os.makedirs(stats_GC_path)
    #Qcycle_stats_path = mapp_stats_path + "/Qcycle"
    #if not os.path.exists(Qcycle_stats_path):
    #    os.makedirs(Qcycle_stats_path)
    #Qdist_stats_path = mapp_stats_path + "/Qdist"
    #if not os.path.exists(Qdist_stats_path):
    #    os.makedirs(Qdist_stats_path)
    #dupMtrx_stats_path = mapp_stats_path + "/DupMatrix"
    #if not os.path.exists(dupMtrx_stats_path):
    #    os.makedirs(dupMtrx_stats_path)
    #summ_stats_path = mapp_stats_path + "/Summary"
    #if not os.path.exists(summ_stats_path):
    #    os.makedirs(summ_stats_path)
    #fig_path = "/mnt/Cromosoma/lmorales/Public/Ustilago/B1/analysis/figures"
    #fig_GCpdf_path = fig_path + "/GCBias_pdf"
    #fig_Qcycle_path = fig_path + "/Qcycle_pdf"
    #fig_Qdist_path = fig_path + "/Qdist_pdf"
    #if not os.path.exists(fig_path):
    #    os.makedirs(fig_path)
    #if not os.path.exists(fig_GCpdf_path):
    #    os.makedirs(fig_GCpdf_path)
    #if not os.path.exists(fig_Qcycle_path):
    #    os.makedirs(fig_Qcycle_path)
    #if not os.path.exists(fig_Qdist_path):
    #    os.makedirs(fig_Qdist_path)
    # Names
    #name_1 = lng_name +"_"+read_1+".good.fq"
    #name_2 = lng_name +"_"+read_2+".good.fq"
    #clean_R1 = fastq_path + "/" + name_1
    #clean_R2 = fastq_path + "/" + name_2
    #bam_name = shrt_name + "_" + map_tool + ".bam"
    #bam_stat_name = shrt_name + "_" + map_tool + "_bam_0_status.txt"
    #bam_stat1_name = shrt_name + "_" + map_tool + "_bam_1_status.txt"
    #bam_stat2_name = shrt_name + "_" + map_tool + "_bam_2_status.txt"
    #bam_mrkdup_name = shrt_name + "_" + map_tool + ".mrkdup.bam"
    #bam_addgp_name = shrt_name + "_" + map_tool + ".mrkdup.addgp.bam"
    #GCBias_name = shrt_name + "_" + map_tool + "_GCBias.txt"
    #GCBias_pdf = shrt_name + "_" + map_tool + "_GCBias.pdf"
    #smmry_name = shrt_name + "_" + map_tool + "_summary_metrics.txt"
    #Qcyc_name = shrt_name + "_" + map_tool + "_Qcycle.txt"
    #Qcyc_pdf = shrt_name + "_" + map_tool + "_Qcycle.pdf"
    #Qdist_name = shrt_name + "_" + map_tool + "_Qdist.txt"
    #Qdist_pdf = shrt_name + "_" + map_tool + "_Qdist.pdf"
    #dupMtrx = shrt_name + "_" + map_tool + "_duplicateMatrix"
    ## Insert size requeriments
    # paths
    #fig_InsSz_pdf_path = "/mnt/Cromosoma/lmorales/Public/Ustilago/B1/analysis/figures/InsertSize_pdf"
    #if not os.path.exists(fig_InsSz_pdf_path):
    #    os.makedirs(fig_InsSz_pdf_path)
    #fig_InsSz_png_path = "/mnt/Cromosoma/lmorales/Public/Ustilago/B1/analysis/figures/InsertSize_png"
    #if not os.path.exists(fig_InsSz_png_path):
    #    os.makedirs(fig_InsSz_png_path)
    #stats_InsSz_path = stats_path + "/mapp/InsertSize_metrics"
    #if not os.path.exists(stats_InsSz_path):
    #    os.makedirs(stats_InsSz_path)
    #r_log_path = "/mnt/Cromosoma/lmorales/Public/Ustilago/B1/log/R/insert_size"
    #if not os.path.exists(r_log_path):
    #    os.makedirs(r_log_path)
    #R_script_path = "/mnt/Cromosoma/lmorales/Public/ymez/bin/scripts/03_mapping"
    # names
    #insrt_name = shrt_name + "_insert_metrics.txt"
    #hist_pdf = shrt_name + "_insert_histogram.pdf"
    #r_log_name = shrt_name + "_histogram.Rout"
    #hist_png = shrt_name + "_insert_histogram.png"
    #print ('''start=$(date +%s.%N)''', file = sge)
    #print ("## MAPPING ", file = sge )
    #if map_tool == "BWA" or map_tool == "bwa":
    #    print ("bwa mem -M -t10 " + reference_genome + " " + clean_R1 + " " + clean_R2 + " | samtools view -hbS - | samtools sort -o " + bam_path + "/" + bam_name + " -", file = sge)
    #elif map_tool == "B2" or map_tool == "Bowtie2":
    #    print ("bowtie2 -t --no-mixed --threads 10 -x " + reference_genome + " -1 " + clean_R1 + " -2 " + clean_R2 + " -S " + sam_path + "/" + sam_name, file = sge)
    #else:
    #    sys.exit("No mapping tool selected, bye-bye")
    print ("#", file = sge)
    print ("this program works in the next directory:" + directory, file = sge)
    print ("the reference is located in:" + wd_ref, file = sge)
    # print ("picard ValidateSamFile I=" + bam_path + "/" + bam_name + " MODE=SUMMARY O=" + bam_stats_path + "/" + bam_stat_name, file = sge)
    # print ("#", file = sge)
    # print ("picard CollectGcBiasMetrics R=" + reference_genome + " I=" + bam_path + "/" + bam_name + " O=" + stats_GC_path + "/" + GCBias_name + " CHART=" + fig_GCpdf_path + "/" + GCBias_pdf + " ASSUME_SORTED=true SUMMARY_OUTPUT=" + summ_stats_path + "/" + smmry_name + " VALIDATION_STRINGENCY=LENIENT", file = sge)
    # print ("#", file = sge)
    # print ("picard MeanQualityByCycle R=" + reference_genome + " I=" + bam_path + "/" + bam_name + " O=" + Qcycle_stats_path + "/" + Qcyc_name + " CHART=" + fig_Qcycle_path + "/" + Qcyc_pdf + " VALIDATION_STRINGENCY=LENIENT", file = sge)
    # print ("#", file = sge)
    # print ("picard QualityScoreDistribution R=" + reference_genome + " I=" + bam_path + "/" + bam_name + " O=" + Qdist_stats_path + "/" + Qdist_name + " CHART=" + fig_Qdist_path + "/" + Qdist_pdf + " VALIDATION_STRINGENCY=LENIENT", file =sge)
    # print ("#", file = sge)
    # print ("picard MarkDuplicates INPUT=" + bam_path + "/" + bam_name + " OUTPUT=" + bam_path + "/" + bam_mrkdup_name + " METRICS_FILE=" + dupMtrx_stats_path + "/" + dupMtrx + " VALIDATION_STRINGENCY=LENIENT", file = sge)
    # print ("#", file = sge)
    # print ("picard ValidateSamFile I=" + bam_path + "/" + bam_mrkdup_name + " MODE=SUMMARY O=" + bam_stats_path + "/" + bam_stat1_name, file = sge)
    # print ("#", file = sge)
    # print ("picard AddOrReplaceReadGroups I=" + bam_path + "/" + bam_mrkdup_name + " O=" + bam_path + "/" + bam_addgp_name + " LB=" + shrt_name + " PL=illumina PU=" + shrt_name + " SM=" + shrt_name + " VALIDATION_STRINGENCY=LENIENT", file = sge)
    # print ("#", file = sge)
    # print ("picard ValidateSamFile I=" + bam_path + "/" + bam_addgp_name + " MODE=SUMMARY O=" + bam_stats_path + "/" + bam_stat2_name, file = sge)
    # print ("#", file = sge)
    # print ("samtools index " + bam_path + "/" + bam_addgp_name, file = sge)
    # print ("#", file = sge)
    # print ("if [[ -s " + bam_path + "/" + bam_addgp_name + " ]]; then", file = sge)
    # print ("rm -f " + bam_path + "/" + bam_mrkdup_name, file = sge)
    # print ("fi", file = sge)
    # print ("## END_MAPPING", file = sge)
    # print ("#", file = sge)
    # print ("#", file = sge)
    # print ("## INSERT SIZE", file = sge )
    # print ("#", file = sge)
    # print ("picard CollectInsertSizeMetrics I=" + bam_path + "/" + bam_addgp_name + " O=" + stats_InsSz_path + "/" + insrt_name + " H=" + fig_InsSz_pdf_path + "/" + hist_pdf + " M=0.5", file = sge)
    # print ("#", file = sge)
    # print ('''R CMD BATCH --no-save --no-restore "--args FILE=''' + "'" + stats_InsSz_path + "/" + insrt_name + "'" + " FILE_OUT=" + "'" + fig_InsSz_png_path + "/" + hist_png + "' bamName='" + bam_path + "/" + bam_addgp_name + "'" + '''" ''' + R_script_path + "/insert_histogram.R " + r_log_path + "/" + r_log_name, file = sge)
    # print ("#", file = sge)
    # print ("## END INSERT SIZE", file = sge)
    # print ('''duration=$(echo "$(date +%s.%N) - $start" | bc)''', file = sge)
    # print ('''execution_time=`printf "%.2f seconds" $duration`''', file = sge)
    # print ('''echo "Script Execution Time: $execution_time"''', file = sge)
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
pyout_name = wd_project + "/log/mapping/"+fecha+"/python/" + tiempo + "_Mapping_py.out"
if not os.path.exists(wd_project + "/log/mapping/"+fecha+"/python/"):
            os.makedirs(wd_project + "/log/mapping/"+fecha+"/python/")
pyout_name = os.path.normpath(pyout_name)
pyoutput = open(pyout_name, "a+")
for i in range(0, len(smpls_ID)):
    ID = smpls_ID[i]    
    #read_1 = extcol(matrix_csv, "Read1")[i]
    #read_2 = extcol(matrix_csv, "Read2")[i]
    #ref_name = extcol(reference_csv, "RefName")[0]
    sge = pipeline(ID)
    #print(tiempo, file = pyoutput)
    #subprocess.run(["qsub",sge], stdout=pyoutput, stderr=subprocess.STDOUT, shell=False, cwd=None, timeout=None, check=True, encoding=None, errors=None, text=None, env=None, universal_newlines=None)
    #print(tiempo, file = pyoutput)
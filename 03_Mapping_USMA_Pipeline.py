###
###
###############################################################################################################################################
### File:           03_Mapping_USMA_Pipeline.py
### Written by:     Jorge Luis Cuamatzi Flores
### Date:           2020_May_29
### Update:         2021_August_13
###
### Project:        This script writes a sge file to map and do the variant calling
### Input:          A csv file with information and simple paths
### Output:         SGE to do mapping and variant calling for several samples
###
# How execute this script: python3 /mnt/Timina/lmorales/Public/Ustilago/C1/bin/scripts/03_Mapping_USMA_Pipeline.py -d /mnt/Timina/lmorales/Public/Ustilago/C1/ -f /mnt/Timina/lmorales/Public/Ustilago/C1/ID.csv -t 03_Mapping -r /mnt/Timina/lmorales/Public/Ustilago/reference/USMA_521_v2.csv -M 8 -w 200
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
ag.add_argument("-d", "--directory", default = "/home/jcuamatzi", help = "path to the project directory")
ag.add_argument("-t", "--task", default = "", help = "Task")
ag.add_argument("-M", "--memory", default = "2", help = "RAM memory")
ag.add_argument("-w", "--window",  default = "100", help = "size of the windows in bp")
##
#
args = vars(ag.parse_args())
arg_file = args["file"]
arg_reference = args["ref"]
email = str(args["email"])
map_tool = str(args["mapping"])
directory = str(args["directory"])
task = str(args["task"])
memory = str(args["memory"])
window = str(args["window"])
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
#$ -pe openmp 10
#$ -m e
source /etc/bashrc
## notification
#$ -M ''' + email + '''
##
## Modules''', file = sge)
    print ("module load bwa/0.7.4 htslib/1.2.1 gcc/5.1.0 samtools/1.9 picard/2.6.0 r/3.6.1", file = sge)
    print ("##", file = sge)
    print ("##", file = sge)

def mapping (smpls_ID):
    save_sge = wd_project + "bin/SGE/" + task + "/"
    if not os.path.exists(save_sge):    
        os.makedirs(save_sge)
    sge_name = save_sge + "/" + fecha + "_" + smpls_ID + "_" + sample_name + ".sge"
    sge = open(sge_name, "w")
    header(smpls_ID,sge)
    reference_genome = wd_ref + "/" + ref_name + ".fa"
    ## Mapping requirements
    ## Paths
    fastq_clean_path = wd_project + "data/fastq/clean/"
    bam_path = wd_project + "data/bam/"
    bam_stats_path = bam_path + "stats/"
    if not os.path.exists(bam_path):
        os.makedirs(bam_path)
    if not os.path.exists(bam_stats_path):
        os.makedirs(bam_stats_path)
    stats_path = wd_project + "analysis/stats/"
    mapp_stats_path = stats_path + "mapp/"
    if not os.path.exists(mapp_stats_path):
        os.makedirs(mapp_stats_path)
    stats_GC_path = mapp_stats_path + "GCBias/"
    if not os.path.exists(stats_GC_path):
        os.makedirs(stats_GC_path)
    Qcycle_stats_path = mapp_stats_path + "Qcycle/"
    if not os.path.exists(Qcycle_stats_path):
        os.makedirs(Qcycle_stats_path)
    Qdist_stats_path = mapp_stats_path + "Qdist/"
    if not os.path.exists(Qdist_stats_path):
        os.makedirs(Qdist_stats_path)
    dupMtrx_stats_path = mapp_stats_path + "DupMatrix/"
    if not os.path.exists(dupMtrx_stats_path):
        os.makedirs(dupMtrx_stats_path)
    summ_stats_path = mapp_stats_path + "Summary/"
    if not os.path.exists(summ_stats_path):
        os.makedirs(summ_stats_path)
    fig_path = wd_project + "analysis/figures/"
    fig_GCpdf_path = fig_path + "GCBias_pdf/"
    fig_Qcycle_path = fig_path + "Qcycle_pdf/"
    fig_Qdist_path = fig_path + "Qdist_pdf/"
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    if not os.path.exists(fig_GCpdf_path):
        os.makedirs(fig_GCpdf_path)
    if not os.path.exists(fig_Qcycle_path):
        os.makedirs(fig_Qcycle_path)
    if not os.path.exists(fig_Qdist_path):
        os.makedirs(fig_Qdist_path)
    # Names
    fastq_name_clean_R1 = smpls_ID + "_" + sample_name + "_R1_clean.fastq.gz"
    fastq_name_clean_R2 = smpls_ID + "_" + sample_name + "_R2_clean.fastq.gz"
    bam_name = smpls_ID + "_" + sample_name + "_" + map_tool + ".bam"
    bam_stat_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_bam_0_status.txt"
    bam_stat1_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_bam_1_status.txt"
    bam_stat2_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_bam_2_status.txt"
    bam_mrkdup_name = smpls_ID + "_" + sample_name + "_" + map_tool + ".mrkdup.bam"
    bam_addgp_name = smpls_ID + "_" + sample_name + "_" + map_tool + ".mrkdup.addgp.bam"
    GCBias_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_GCBias.txt"
    GCBias_pdf = smpls_ID + "_" + sample_name + "_" + map_tool + "_GCBias.pdf"
    smmry_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_summary_metrics.txt"
    Qcyc_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_Qcycle.txt"
    Qcyc_pdf = smpls_ID + "_" + sample_name + "_" + map_tool + "_Qcycle.pdf"
    Qdist_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_Qdist.txt"
    Qdist_pdf = smpls_ID + "_" + sample_name + "_" + map_tool + "_Qdist.pdf"
    dupMtrx = smpls_ID + "_" + sample_name + "_" + map_tool + "_duplicateMatrix"
    #
    #Files
    clean_R1 = fastq_clean_path + fastq_name_clean_R1
    clean_R2 = fastq_clean_path + fastq_name_clean_R2
    #
    ## Insert size requeriments
    # Insert size requeriments: paths
    fig_InsSz_pdf_path = fig_path + "InsertSize_pdf/"
    if not os.path.exists(fig_InsSz_pdf_path):
        os.makedirs(fig_InsSz_pdf_path)
    fig_InsSz_png_path = fig_path + "InsertSize_png/"
    if not os.path.exists(fig_InsSz_png_path):
        os.makedirs(fig_InsSz_png_path)
    stats_InsSz_path = stats_path + "mapp/InsertSize_metrics/"
    if not os.path.exists(stats_InsSz_path):
        os.makedirs(stats_InsSz_path)
    r_log_path = wd_project + "log/R/insert_size/"
    if not os.path.exists(r_log_path):
        os.makedirs(r_log_path)
    R_script_path = "/mnt/Timina/lmorales/Public/ymez/bin/scripts/03_mapping"
    # Insert size requeriments: names
    insrt_name = smpls_ID + "_insert_metrics.txt"
    hist_pdf = smpls_ID + "_insert_histogram.pdf"
    r_log_name = smpls_ID + "_histogram.Rout"
    hist_png = smpls_ID + "_insert_histogram.png"
    ## Coverage requirements
    # Coverage requirements: paths
    cov_path = wd_project + "analysis/Coverage/"
    if not os.path.exists(cov_path):
        os.makedirs(cov_path)
    # Coverage requirements: names
    cov_name = smpls_ID + "_" + sample_name + ".coverage_per.bp"
    cov_name_bed5 = smpls_ID + "_" + sample_name + ".coverage.bed5"
    cov_name_wind = smpls_ID + "_" + sample_name + ".coverage_" + window + "bp_windows.bed"
    #
    print ('''start=$(date +%s.%N)''', file = sge)
    print ("## MAPPING ", file = sge )
    print ("bwa mem -M -t10 " + reference_genome + " " + clean_R1 + " " + clean_R2 + " | samtools view -hbS - | samtools sort -o " + bam_path + bam_name + " -", file = sge)
    print ("#", file = sge)
    print ("picard ValidateSamFile I=" + bam_path + bam_name + " MODE=SUMMARY O=" + bam_stats_path + bam_stat_name, file = sge)
    print ("#", file = sge)
    print ("picard CollectGcBiasMetrics R=" + reference_genome + " I=" + bam_path + bam_name + " O=" + stats_GC_path + GCBias_name + " CHART=" + fig_GCpdf_path + GCBias_pdf + " ASSUME_SORTED=true SUMMARY_OUTPUT=" + summ_stats_path + "/" + smmry_name + " VALIDATION_STRINGENCY=LENIENT", file = sge)
    print ("#", file = sge)
    print ("picard MeanQualityByCycle R=" + reference_genome + " I=" + bam_path + bam_name + " O=" + Qcycle_stats_path + Qcyc_name + " CHART=" + fig_Qcycle_path + Qcyc_pdf + " VALIDATION_STRINGENCY=LENIENT", file = sge)
    print ("#", file = sge)
    print ("picard QualityScoreDistribution R=" + reference_genome + " I=" + bam_path + bam_name + " O=" + Qdist_stats_path + Qdist_name + " CHART=" + fig_Qdist_path + Qdist_pdf + " VALIDATION_STRINGENCY=LENIENT", file =sge)
    print ("#", file = sge)
    print ("picard MarkDuplicates INPUT=" + bam_path + bam_name + " OUTPUT=" + bam_path + bam_mrkdup_name + " METRICS_FILE=" + dupMtrx_stats_path + dupMtrx + " VALIDATION_STRINGENCY=LENIENT", file = sge)
    print ("#", file = sge)
    print ("picard ValidateSamFile I=" + bam_path + bam_mrkdup_name + " MODE=SUMMARY O=" + bam_stats_path + bam_stat1_name, file = sge)
    print ("#", file = sge)
    print ("picard AddOrReplaceReadGroups I=" + bam_path + bam_mrkdup_name + " O=" + bam_path + bam_addgp_name + " LB=" + smpls_ID + " PL=illumina PU=" + smpls_ID + " SM=" + smpls_ID + " VALIDATION_STRINGENCY=LENIENT", file = sge)
    print ("#", file = sge)
    print ("picard ValidateSamFile I=" + bam_path + bam_addgp_name + " MODE=SUMMARY O=" + bam_stats_path + bam_stat2_name, file = sge)
    print ("#", file = sge)
    print ("samtools index " + bam_path + bam_addgp_name, file = sge)
    print ("#", file = sge)
    print ("if [[ -s " + bam_path + bam_addgp_name + " ]]; then", file = sge)
    print ("rm -f " + bam_path + bam_mrkdup_name, file = sge)
    print ("fi", file = sge)
    print ("## END_MAPPING", file = sge)
    print ("#", file = sge)
    print ("#", file = sge)
    ## Insert size
    print ("## INSERT SIZE", file = sge )
    print ("#", file = sge)
    print ("picard CollectInsertSizeMetrics I=" + bam_path + bam_addgp_name + " O=" + stats_InsSz_path + insrt_name + " H=" + fig_InsSz_pdf_path + hist_pdf + " M=0.5", file = sge)
    print ("#", file = sge)
    print ('''R CMD BATCH --no-save --no-restore "--args FILE=''' + "'" + stats_InsSz_path + insrt_name + "'" + " FILE_OUT=" + "'" + fig_InsSz_png_path + hist_png + "' bamName='" + bam_path + bam_addgp_name + "'" + '''" ''' + R_script_path + "/insert_histogram.R " + r_log_path + r_log_name, file = sge)
    print ("#", file = sge)
    print ("## END INSERT SIZE", file = sge)
    print ("#", file = sge)
    print ("#", file = sge)
    ## Coverage
    print ("## COVERAGE", file = sge)
    print ("## Obtain genome coverage in " + window + " bp non-overlapping windows", file = sge)
    print ("bedtools genomecov -ibam " + bam_path + bam_addgp_name + " -d > " + cov_path + cov_name, file = sge)
    print ('''awk -vFS="\\t" -vOFS="\\t" '{ print $1, $2, ($2 +1), ".", $3 }' ''' + cov_path + cov_name + " | sort-bed - > " + cov_path + cov_name_bed5, file = sge)
    print ("bedops --merge " + cov_path + cov_name_bed5 + " | bedops --chop " + window + " - | bedmap --echo --mean --delim '\\t' - " + cov_path + cov_name_bed5 + " > " + cov_path + cov_name_wind, file = sge)
    print ("## END COVERAGE", file = sge)
    print ('''duration=$(echo "$(date +%s.%N) - $start" | bc)''', file = sge)
    print ('''execution_time=`printf "%.2f seconds" $duration`''', file = sge)
    print ('''echo "Script Execution Time: $execution_time"''', file = sge)
    sge.close()
    return sge_name  
# Update 16/08/2021

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
    #print(tiempo, file = pyoutput)
    #subprocess.run(["qsub",sge], stdout=pyoutput, stderr=subprocess.STDOUT, shell=False, cwd=None, timeout=None, check=True, encoding=None, errors=None, text=None, env=None, universal_newlines=None)
    #print(tiempo, file = pyoutput)
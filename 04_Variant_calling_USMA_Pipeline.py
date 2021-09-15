###
###
###############################################################################################################################################
### File:           04_Variant_calling_USMA_Pipeline.py
### Written by:     Jorge Luis Cuamatzi Flores
### Date:           2021_August_16
###
### Project:        This script writes a sge file to map and do the variant calling
### Input:          A csv file with information and simple paths
### Output:         SGE to do mapping and variant calling for several samples
###
# How execute this script: python3 /mnt/Timina/lmorales/Public/Ustilago/C1/bin/scripts/04_Variant_calling_USMA_Pipeline.py -d /mnt/Timina/lmorales/Public/Ustilago/C1/ -f /mnt/Timina/lmorales/Public/Ustilago/C1/ID.csv -t 04_Variant_calling -r /mnt/Timina/lmorales/Public/Ustilago/reference/USMA_521_v2.csv -M 16 -w 200
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
ag.add_argument("-g", "--gff",  default = "", help = "size of the windows in bp")
##
#
args = vars(ag.parse_args())
arg_file = args["file"]
reference_genome = str(args["ref"])
email = str(args["email"])
map_tool = str(args["mapping"])
directory = str(args["directory"])
task = str(args["task"])
memory = str(args["memory"])
window = str(args["window"])
gff_file = str(args["gff"])
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
    print ("module load gcc/5.1.0 bedops/2.4.20 bedtools/2.24 gatk/4.1.1.0 picard/2.6.0 samtools/1.9 bcftools/1.9 vcftools/0.1.14" , file = sge)
    print ("##", file = sge)
    print ("##", file = sge)

def variantcall (smpls_ID):
    save_sge = wd_project + "bin/SGE/" + task + "/"
    if not os.path.exists(save_sge):    
        os.makedirs(save_sge)
    sge_name = save_sge + "/" + fecha + "_" + smpls_ID + "_" + sample_name + ".sge"
    sge = open(sge_name, "w")
    header(smpls_ID,sge)
    ## Paths
    bam_path = wd_project + "data/bam/"
    vcf_path = wd_project + "data/VCF/"
    if not os.path.exists(vcf_path):
        os.makedirs(vcf_path)
    g_vcf_path = vcf_path + "g_vcf/"
    if not os.path.exists(g_vcf_path):
        os.makedirs(g_vcf_path)
    gt_vcf_path = vcf_path + "gt_vcf/"
    if not os.path.exists(gt_vcf_path):
        os.makedirs(gt_vcf_path)
    SNP_vcf_path = gt_vcf_path + "SNP/"
    if not os.path.exists(SNP_vcf_path):
        os.makedirs(SNP_vcf_path)
    SNP_flt_vcf_path = SNP_vcf_path + "filter_vcf/"
    if not os.path.exists(SNP_flt_vcf_path):
        os.makedirs(SNP_flt_vcf_path)
    SNP_flt_vcf_PASS_path = SNP_flt_vcf_path + "PASS/"
    if not os.path.exists(SNP_flt_vcf_PASS_path):
        os.makedirs(SNP_flt_vcf_PASS_path)
    INDEL_vcf_path = gt_vcf_path + "INDEL/"
    if not os.path.exists(INDEL_vcf_path):
        os.makedirs(INDEL_vcf_path)    
    INDEL_flt_vcf_path = INDEL_vcf_path + "filter_vcf/"
    if not os.path.exists(INDEL_flt_vcf_path):
        os.makedirs(INDEL_flt_vcf_path)
    INDEL_flt_vcf_PASS_path = INDEL_flt_vcf_path + "PASS/"
    if not os.path.exists(INDEL_flt_vcf_PASS_path):
        os.makedirs(INDEL_flt_vcf_PASS_path)
    
    annotation_vcf_path = vcf_path + "Annotation/"
    if not os.path.exists(annotation_vcf_path):
        os.makedirs(annotation_vcf_path)
    annttd_SNP_path = annotation_vcf_path + "SNP/"
    if not os.path.exists(annttd_SNP_path):
        os.makedirs(annttd_SNP_path)
    annttd_BG_SNP_path = annotation_vcf_path + "Background/SNP/"
    if not os.path.exists(annttd_BG_SNP_path):
        os.makedirs(annttd_BG_SNP_path)
    annttd_BG_INDEL_path = annotation_vcf_path + "Background/INDEL/"
    if not os.path.exists(annttd_BG_INDEL_path):
        os.makedirs(annttd_BG_INDEL_path)    
    annttd_SNP_raw_table_path = annttd_SNP_path + "Tables/SNP_raw/"
    if not os.path.exists(annttd_SNP_raw_table_path):
        os.makedirs(annttd_SNP_raw_table_path)
    annttd_SNP_flt_table_path = annttd_SNP_path + "Tables/SNP_filter/"
    if not os.path.exists(annttd_SNP_flt_table_path):
        os.makedirs(annttd_SNP_flt_table_path)
    annttd_SNP_PASS_table_path = annttd_SNP_path + "Tables/SNP_PASS/"
    if not os.path.exists(annttd_SNP_PASS_table_path):
        os.makedirs(annttd_SNP_PASS_table_path)
    annttd_SNP_BG_table_path = annttd_SNP_path + "Tables/SNP_Background/"
    if not os.path.exists(annttd_SNP_BG_table_path):
        os.makedirs(annttd_SNP_BG_table_path)
    annttd_INDEL_path = annotation_vcf_path + "INDEL/"
    if not os.path.exists(annttd_INDEL_path):
        os.makedirs(annttd_INDEL_path)
    annttd_INDELS_BG_table_path = annttd_INDEL_path + "Tables/INDEL_Background/"
    if not os.path.exists(annttd_INDELS_BG_table_path):
        os.makedirs(annttd_INDELS_BG_table_path)
    annttd_INDEL_raw_table_path = annttd_INDEL_path + "Tables/INDEL_raw/"
    if not os.path.exists(annttd_INDEL_raw_table_path):
        os.makedirs(annttd_INDEL_raw_table_path)
    annttd_INDEL_flt_table_path = annttd_INDEL_path + "Tables/INDEL_filter/"
    if not os.path.exists(annttd_INDEL_flt_table_path):
        os.makedirs(annttd_INDEL_flt_table_path)
    annttd_INDEL_PASS_table_path = annttd_INDEL_path + "Tables/INDEL_PASS/"
    if not os.path.exists(annttd_INDEL_PASS_table_path):
        os.makedirs(annttd_INDEL_PASS_table_path)
    
    
    
    
    
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
    
    
    bam_name = smpls_ID + "_" + sample_name + "_" + map_tool + ".mrkdup.addgp.bam"
    bam_name_2 = smpls_ID + "_" + sample_name + "_bamout_VC.bam"
    g_vcf_name = smpls_ID + "_" + sample_name + "_" + map_tool + ".g.vcf.gz" # compress
    gt_vcf_name = smpls_ID + "_" + sample_name + "_" + map_tool + ".gt.g.vcf.gz"
    SNP_vcf_name = smpls_ID + "_" + sample_name + "_" + map_tool + ".gt.SNP.g.vcf.gz"
    SNP_flt_vcf_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_SNP_flt.vcf.gz"
    SNP_flt_vcf_PASS_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_SNP_flt_PASS.vcf.gz"
    INDEL_vcf_name = smpls_ID + "_" + sample_name + "_" + map_tool + ".gt.INDEL.g.vcf.gz"
    INDEL_flt_vcf_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_INDEL_flt.vcf.gz"
    INDEL_flt_vcf_PASS_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_INDEL_flt_PASS.vcf.gz"
    SNP_flt_vcf_gz_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_SNP_flt.vcf.gz"
    INDEL_flt_vcf_gz_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_INDEL_flt.vcf.gz"
    SNP_flt_vcf_gz_PASS_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_SNP_flt_PASS.vcf.gz"
    INDEL_flt_vcf_gz_PASS_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_INDEL_flt_PASS.vcf.gz"
    SNP_flt_vcf_PASS_BG_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_SNP_PASS_BG.vcf.gz"
    INDEL_flt_vcf_PASS_BG_name = smpls_ID + "_" + sample_name + "_" + map_tool + "_INDEL_PASS_BG.vcf.gz"
    
    annttd_SNP_name = smpls_ID + "_" + sample_name + "_SNP_raw_annotated.vcf.gz"
    annttd_INDEL_name = smpls_ID + "_" + sample_name + "_INDEL_raw_annotated.vcf.gz"
    annttd_SNP_flt_name = smpls_ID + "_" + sample_name + "_SNP_flt_annotated.vcf.gz"
    annttd_INDEL_flt_name = smpls_ID + "_" + sample_name + "_INDEL_flt_annotated.vcf.gz"
    annttd_SNP_PASS_name = smpls_ID + "_" + sample_name + "_SNP_PASS_annotated.vcf.gz"
    annttd_INDEL_PASS_name = smpls_ID + "_" + sample_name + "_INDEL_PASS_annotated.vcf.gz"
    annttd_SNP_PASS_BG_name = smpls_ID + "_" + sample_name + "_SNP_PASS_BG_annotated.vcf.gz"
    annttd_INDEL_PASS_BG_name = smpls_ID + "_" + sample_name + "_INDEL_PASS_BG_annotated.vcf.gz"
    
    bamout = bam_path + "bamout_VC/"
    if not os.path.exists(bamout):
        os.makedirs(bamout)
    
    #
    print ('''start=$(date +%s.%N)''', file = sge)
    print ("## VARIANT CALLING", file = sge )
    print ("#", file = sge)
    print ('''gatk --java-options "-Xmx16g" HaplotypeCaller -R ''' + reference_genome + " -I " + bam_path + bam_addgp_name + " -O " + g_vcf_path + g_vcf_name + " -bamout " + bamout + bam_name_2 + " --sample-ploidy 1 --annotation DepthPerSampleHC --annotation StrandBiasBySample --annotation AlleleFraction --annotation AS_FisherStrand --annotation ChromosomeCounts --dont-use-soft-clipped-bases TRUE -ERC GVCF", file = sge)
    print ("#", file = sge)
    print ('''gatk --java-options "-Xmx16g" GenotypeGVCFs -R ''' + reference_genome + " -V " + g_vcf_path + g_vcf_name + " -O " + gt_vcf_path + gt_vcf_name + " --sample-ploidy 1 --annotation DepthPerSampleHC --annotation StrandBiasBySample --annotation AlleleFraction --annotation AS_FisherStrand --annotation ChromosomeCounts", file = sge)
    print ("#", file = sge)
    print ('''gatk --java-options "-Xmx16g" SelectVariants -R ''' + reference_genome + " -V " + gt_vcf_path + gt_vcf_name + " --select-type-to-include SNP -O " + SNP_vcf_path + SNP_vcf_name, file = sge)
    print ("#", file = sge)
    print ("gatk VariantFiltration -R " + reference_genome + " -V " + SNP_vcf_path + SNP_vcf_name + " -O " + SNP_flt_vcf_path + SNP_flt_vcf_name + ''' --filter-name "SNP_QD_FS_MQ_SOR_filters" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0"''', file = sge)
    print ("#", file = sge)
    print ('''awk '/^#/||$7=="PASS"' ''' + SNP_flt_vcf_path + SNP_flt_vcf_name + " > " + SNP_flt_vcf_PASS_path + SNP_flt_vcf_PASS_name, file = sge)
    print ("#", file = sge)
    print ('''gatk --java-options "-Xmx16g" SelectVariants -R ''' + reference_genome + " -V " + gt_vcf_path + gt_vcf_name + " --select-type-to-include INDEL -O " + INDEL_vcf_path + INDEL_vcf_name, file = sge)
    print ("#", file = sge)
    print ("gatk VariantFiltration -R " + reference_genome + " -V " + INDEL_vcf_path + INDEL_vcf_name + " -O " + INDEL_flt_vcf_path + INDEL_flt_vcf_name + ''' --filter-name "INDEL_QD_FS_Q_filters" --filter-expression "QD < 2.0 || FS > 200.0 || QUAL < 30.0"''', file = sge)
    print ("#", file = sge)
    print ('''awk '/^#/||$7=="PASS"' ''' + INDEL_flt_vcf_path + INDEL_flt_vcf_name + " > " + INDEL_flt_vcf_PASS_path + INDEL_flt_vcf_PASS_name, file = sge)
    print ("#", file = sge)
    print ("## Annotating the files", file = sge)
    print ("# raw vcf files", file = sge)
    print ("# SNP", file = sge)
    print ("bcftools csq -f " + reference_genome + " -g " + gff_file + SNP_vcf_path + SNP_vcf_name + " -Oz -o " + annttd_SNP_path + annttd_SNP_name, file = sge)
    print ("bcftools query -f'[%CHROM\\t%POS\\t%SAMPLE\\t%TBCSQ\\n]' " + annttd_SNP_path + annttd_SNP_name + " > " + annttd_SNP_raw_table_path + smpls_ID + "_" + sample_name + "_raw_AnnTable.txt", file = sge)
    print ('''awk -F"," '{print $1"\\t"$3}' ''' + annttd_SNP_raw_table_path + smpls_ID + "_" + sample_name + '''_raw_AnnTable.txt | awk -F"|" '{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7}' > ''' + annttd_SNP_raw_table_path + smpls_ID + "_" + sample_name + "_raw_AnnTable.txt2", file = sge)
    print ("rm " + annttd_SNP_raw_table_path + smpls_ID + "_" + sample_name + "_raw_AnnTable.txt", file = sge)
    print ("cp " + annttd_SNP_raw_table_path + smpls_ID + "_" + sample_name + "_raw_AnnTable.txt2 " + annttd_SNP_raw_table_path + smpls_ID + "_" + sample_name + "_raw_AnnTable.txt", file = sge)
    print ("rm " + annttd_SNP_raw_table_path + smpls_ID + "_" + sample_name + "_raw_AnnTable.txt2", file = sge)
    print ("#", file = sge) 
    print ("# INDEL", file = sge)
    print ("bcftools csq -f " + genom_ref + " -g " + gff_path + "/USMA_521_v2_allChr.gff " + INDEL_vcf_path + "/" + INDEL_vcf_name + ".gz -Oz -o " + annttd_INDEL_path + "/" + annttd_INDEL_name, file = sge)   
    print ("bcftools query -f'[%CHROM\\t%POS\\t%SAMPLE\\t%TBCSQ\\n]' " + annttd_INDEL_path + "/" + annttd_INDEL_name + " > " + annttd_INDEL_raw_table_path + "/" + shrt_name + "_raw_AnnTable.txt", file = sge)
    print ('''awk -F"," '{print $1"\\t"$2"\\t"$3"\\t"$5"\\t"}' ''' + annttd_INDEL_raw_table_path + "/" + shrt_name + '''_raw_AnnTable.txt | awk -F"|" '{print $1"\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}' > ''' + annttd_INDEL_raw_table_path + "/" + shrt_name + "_raw_AnnTable.txt2", file = sge)
    print ("rm " + annttd_INDEL_raw_table_path + "/" + shrt_name + "_raw_AnnTable.txt", file = sge)
    print ("cp " + annttd_INDEL_raw_table_path + "/" + shrt_name + "_raw_AnnTable.txt2 " + annttd_INDEL_raw_table_path + "/" + shrt_name + "_raw_AnnTable.txt", file = sge)
    print ("rm " + annttd_INDEL_raw_table_path + "/" + shrt_name + "_raw_AnnTable.txt2", file = sge)
    print ("#", file = sge)
    
    print ("", file =sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    ## Insert size
    print ("", file = sge )
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    ## Coverage
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    print ("", file = sge)
    sge.close()
    return sge_name  
# Update 16/08/2021

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
    sge = variantcall(ID)
    #print(tiempo, file = pyoutput)
    #subprocess.run(["qsub",sge], stdout=pyoutput, stderr=subprocess.STDOUT, shell=False, cwd=None, timeout=None, check=True, encoding=None, errors=None, text=None, env=None, universal_newlines=None)
    #print(tiempo, file = pyoutput)
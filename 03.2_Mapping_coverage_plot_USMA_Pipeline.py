###
###
###############################################################################################################################################
### File:           03.2_Mapping_coverage_plot_USMA_Pipeline.py
### Written by:     Jorge Luis Cuamatzi Flores
### Date:           2020_May_29
### Update:         2021_August_13
###
### Project:        This script writes a Rscript file to map and do the variant calling
### Input:          A csv file with information and simple paths
### Output:         SGE to do mapping and variant calling for several samples
###
# How execute this script: python3 /mnt/Timina/lmorales/Public/Ustilago/C1/bin/scripts/03.2_Mapping_coverage_plot_USMA_Pipeline.py -d /mnt/Timina/lmorales/Public/Ustilago/C1/ -f /mnt/Timina/lmorales/Public/Ustilago/C1/ID.csv -t Mapping -r /mnt/Timina/lmorales/Public/Ustilago/reference/USMA_521_v2.csv -M 8 -w 200
###############################################################################################################################################
## Libraries
import argparse
import os
from datetime import date, datetime
import sys
import subprocess
##
fecha = date.today().strftime('%Y%m%d')
tiempo = datetime.now().strftime("%H.%M.%S")
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
# Function to write the header of Rscript file
#
#def header(smpls_ID, Rscript):
def header(smpls_ID,Rscript):
    print ("#libraries", file = Rscript)
    print ("library(scales)", file = Rscript)
    print ("library(ggplot2)", file = Rscript)
    print ("##", file = Rscript)

def mapping (smpls_ID):
    save_Rscript = wd_project + "bin/Rscript/" + task + "/"
    if not os.path.exists(save_Rscript):
        os.makedirs(save_Rscript)
    Rscript_name = save_Rscript + "/" + fecha + "_" + smpls_ID + ".R"
    Rscript = open(Rscript_name, "w")
    header(smpls_ID,Rscript)
    reference_genome = wd_ref + "/" + ref_name + ".fa"
    ## 
    ## Coverage requirements
    # Coverage requirements: paths
    cov_path = wd_project + "analysis/coverage/"
    # Coverage requirements: names
    cov_name_wind = smpls_ID + "_" + sample_name + ".coverage_" + window + "bp_windows.bed"
    plot_path = wd_project + "analysis/figures/Coverage/"
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    plot_name = smpls_ID + "_depth_coverage_" + window + "bp.png"
    #
    print ("#", file = Rscript)
    print ("rm(list=ls())", file = Rscript)
    print ("s." + smpls_ID + ' <- read.delim("' + cov_path + cov_name_wind + '", header = F)', file = Rscript)
    print ("names(s." + smpls_ID + ")[names(s." + smpls_ID + ') == "V1"] <- "Chr"', file = Rscript)
    print ("names(s." + smpls_ID + ")[names(s." + smpls_ID + ') == "V2"] <- "Start"', file = Rscript)
    print ("names(s." + smpls_ID + ")[names(s." + smpls_ID + ') == "V3"] <- "End"', file = Rscript)
    print ("names(s." + smpls_ID + ")[names(s." + smpls_ID + ') == "V4"] <- "Cov"', file = Rscript)
    print ("#", file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_1"] <- "I"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_2"] <- "II"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_3"] <- "III"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_4"] <- "IV"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_5"] <- "V"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_6"] <- "VI"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_7"] <- "VII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_8"] <- "VIII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_9"] <- "IX"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_10"] <- "X"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_11"] <- "XI"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_12"] <- "XII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_13"] <- "XIII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_14"] <- "XIV"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_15"] <- "XV"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_16"] <- "XVI"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_17"] <- "XVII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_18"] <- "XVIII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_19"] <- "XIX"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_20"] <- "XX"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_21"] <- "XXI"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_22"] <- "XXII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_23"] <- "XXIII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_24"] <- "Mt"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_25"] <- "XV"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_26"] <- "XVI"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_27"] <- "XVII"', file = Rscript)
    print ("s." + smpls_ID + '$Chr[s.' + smpls_ID + '$Chr == "USMA_521_v2_28"] <- "XVIII"', file = Rscript)
    print ("s." + smpls_ID + "$Chr <- factor(s." + smpls_ID + '''$Chr, levels = c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI","XVII","XVIII","XIX","XX","XXI","XXII","XXIII","Mt","XXV","XXVI","XXVII","XXVIII"))''', file = Rscript)
    print ('''scientific <- function(x){
            ifelse(x==0, "0",
                        parse(text=gsub("[+]", "",gsub("e", " %*% 10^", scientific_format()(x)))))
                        } ''', file = Rscript)
    print ("#", file = Rscript)
    print ("plot." + smpls_ID + " <- ggplot(s." + smpls_ID + ''', aes(x = End, y = Cov, colour = Chr)) +
  geom_point() +
  facet_grid(.~ Chr, space = 'free_x', scales = "free_x" ) + 
  theme_bw() +
  #scale_y_log10(labels = scientific)+
  scale_y_continuous(limits = c(0,1000))+
  labs(title = "Sample: ''' + smpls_ID + " (" + window + ''' non-overlapping bp)", x = "Chromosome", y = "Media coverage depth") +
  theme(legend.position = "none",
        panel.spacing.x = grid::unit(0, "cm"),
        panel.border = element_rect(colour = "grey",
                                    size = 0.1),
        panel.ontop = FALSE,
        axis.title.x = element_text(face = "bold",
                                    color = "black", 
                                    size = 22),
        axis.title.y = element_text(face = "bold",
                                    color = "black", 
                                    size = 22),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", 
                                   size = 18),
        axis.ticks.x= element_blank(),
        plot.title = element_text(face = "bold", 
                                  color = "red", 
                                  size = 24, 
                                  hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14,
                                  face = "bold"))''', file = Rscript)
    print('''ggsave("''' + plot_path + plot_name + '''", plot = plot.''' + smpls_ID + ''', width = 20, height = 12, units = "in", dpi = 150)''', file = Rscript)
    print ("rm(list=ls())", file = Rscript)
    print ("#", file = Rscript)
    Rscript.close()
    return Rscript_name  
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
    Rscript = mapping(ID)
    #subprocess.run(["Rscript",Rscript], stdout=pyoutput)
    #print(tiempo, file = pyoutput)
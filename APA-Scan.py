import csv
import time
import sys
from operator import itemgetter
import math
import pandas as pd
from Bio import SeqIO
import re
import bisect
from bisect import bisect_left
from scipy.stats import chisquare
import numpy as np
import methods
import glob, os
import preprocess
import configparser

(speciesFlag, inputFlag, outFlag, pasFlag) = (0, 0, 0, 0)

startTime = time.time()

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

config = configparser.ConfigParser()
config.read('configuration.ini')

input1_dir = config['INPUT_RNAseq']['input1']
input2_dir = config['INPUT_RNAseq']['input2']
if input1_dir[-1] == "/":
	input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
	input2_dir = input2_dir[:-1]
pasSeq_dir1 = config['INPUT_PASseq']['pas1']
pasSeq_dir2 = config['INPUT_PASseq']['pas2']
if pasSeq_dir1[:-1] == "/":
	pasSeq_dir1 = pasSeq_dir1[:-1]
if pasSeq_dir2[:-1] == "/":
	pasSeq_dir2 = pasSeq_dir2[:-1]
output_dir = config['OUTPUT_FOLDER']['output_dir']
if output_dir[-1] != "/":
	output_dir += "/"
extended = config['Extended_3UTR']['extended']
all_events = config['All_events']['All']

os.makedirs(output_dir, exist_ok=True)
inp_annotation = config['ANNOTATION']['annotation']
ref_genome = config['ANNOTATION']['genome']
number_of_cores = config['CORES']['number_of_cores']
g1_name = input1_dir.split("/")[-1]
g2_name = input2_dir.split("/")[-1]
print("RNA-seq input 1 dir:", input1_dir)
print("RNA-seq input 2 dir:", input2_dir)
print("3'-end-seq input 1 dir:", pasSeq_dir1)
print("3'-end-seq input 2 dir:", pasSeq_dir2)
print("Output Dir:", output_dir)
print("Number of Cores:", number_of_cores)
print("Annotation:", inp_annotation, ref_genome, "\n\n")
print("Loading chromosomes...")
df = pd.read_csv(inp_annotation, delimiter='\t')
chr_list = df['chrom'].str.split("_", n = 1, expand = True)
chromosomes = list(set(chr_list[0])-set(['chrM', 'chrUn']))

chromDict = methods.makeChromDict(chromosomes, inp_annotation)

print("Creating read coverage files for RNA-seq data...")

os.chdir(input1_dir)
for sample1 in glob.glob("*.bam"):
    preprocess.SamtoText(input1_dir, sample1, chromosomes, number_of_cores)
os.chdir(input2_dir)
for sample2 in glob.glob("*.bam"):
    preprocess.SamtoText(input2_dir, sample2, chromosomes, number_of_cores)

result_filename = "APA_Scan_"+g1_name+"_Vs_"+g2_name
if pasSeq_dir1 == 'NULL' or pasSeq_dir2=='NULL':
	s1_namelist = list_dirs(input1_dir)
	s2_namelist = list_dirs(input2_dir)

	print("Preparing result using RNA-seq data only")
	methods.Get_Signal_Positions(chromosomes, chromDict, inp_annotation, ref_genome, output_dir, extended)
	if all_events == 'all':
		methods.with_PAS_signal_all(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, result_filename)
	else:
		methods.with_PAS_signal(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, result_filename)
	os.remove(output_dir+"Signal_positions.csv")

else:
	print("Creating read coverage files for 3'-end-seq data...")
	os.chdir(pasSeq_dir1)
	for sample1 in glob.glob("*.bam"):
	    preprocess.SamtoText(pasSeq_dir1, sample1, chromosomes, number_of_cores)
	os.chdir(pasSeq_dir2)
	for sample2 in glob.glob("*.bam"):
	    preprocess.SamtoText(pasSeq_dir2, sample2, chromosomes, number_of_cores)

	p1_namelist = list_dirs(pasSeq_dir1)
	p2_namelist = list_dirs(pasSeq_dir2)
	p1_name = pasSeq_dir1.split("/")[-1]
	p2_name = pasSeq_dir2.split("/")[-1]

	methods.Get_Peak_Positions(chromosomes, chromDict, inp_annotation, pasSeq_dir1, pasSeq_dir2, p1_name, p2_name, output_dir, extended)
	if all_events == 'all':
		methods.with_PA_peaks_all(chromosomes, input1_dir, input2_dir, g1_name, g2_name, output_dir, result_filename)
	else:
		methods.with_PA_peaks(chromosomes, input1_dir, input2_dir, g1_name, g2_name, output_dir, result_filename)
	os.remove(output_dir+"Peak_positions.csv")

print("Total time:", round((time.time() - startTime)/60, 2), "minutes.")
read_file = pd.read_csv(output_dir+result_filename+".csv", delimiter = '\t')
read_file.to_excel (output_dir+result_filename+".xlsx", index = None, header=True)
os.remove(output_dir+result_filename+".csv")

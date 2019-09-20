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
import peakutils
import numpy as np
import methods
import os
import preprocess

(speciesFlag, inputFlag, outFlag, pasFlag) = (0, 0, 0, 0)

startTime = time.time()

if(len(sys.argv)<5):
	print("Please provide all of the mandatory arguments. Example: $ python3 APA-scan.py annotation ref_genome input_directory1 input_directory2")
	sys.exit()

# python3 APA-scan.py annotation ref_genome input_directory1 input_directory2 -o output_directory -p pas_seq_dir1 pas_seq_dir2

if len(sys.argv)>=5:
	for ii in range(len(sys.argv)):
		if sys.argv[ii] == '-o' or sys.argv[ii] == '-O':
			outFlag = 1
			output_dir = sys.argv[ii+1]+'/'
		elif sys.argv[ii] == '-p' or sys.argv[ii] == '-P':
			pasFlag = 1
			pasSeq_input1 = sys.argv[ii+1]
			pasSeq_input2 = sys.argv[ii+2]

inp_annotation = sys.argv[1]
ref_genome = sys.argv[2]
input1_dir = sys.argv[3]
input2_dir = sys.argv[4]

if outFlag == 0:
	output_dir = 'Output-APA-Scan/'

os.makedirs(output_dir, exist_ok=True)

s1_dir, s2_dir, samplenames, bamfile_names = methods.processInputFiles(input1_dir, input2_dir)

if len(bamfile_names) == 0:
	print("Bamfiles not found")
	sys.exit()

chromosomes = []
with open(inp_annotation, 'r') as f:
    reader = csv.reader(f, dialect='excel', delimiter='\t')
    headers = next(f)
    annotList = list(reader)
    for rows in annotList:
    	if '_' not in rows[1] and rows[1]!='chrM':
	    	chromosomes.append(rows[1])
    chr_set = set(chromosomes)
    chromosomes = list(chr_set)

#print(chromosomes)

if pasFlag == 0:
	for i in range(2):
		if i==0:
			inputPath = s1_dir
		else:
			inputPath = s2_dir

	preprocess.SamtoText(input_path, samplenames[i], bamfile_names[i], chromosomes)
	filename, result_filename = methods.Generate_withPasSeqSignal(chromosomes, inp_annotation, ref_genome, output_dir)
else:
	p1_dir, p2_dir, pas_samplenames, pas_bamfile_names = methods.processInputFiles(pasSeq_input1, pasSeq_input2)

	if len(pas_bamfile_names) == 0:
		print("Bamfiles not found")
		sys.exit()

	for i in range(2):
		if i==0:
			inputPath = p1_dir
		else:
			inputPath = p2_dir
		preprocess.SamtoText(input_path, pas_samplenames[i], pas_bamfile_names[i], chromosomes)

	filename, result_filename = methods.Generate_withPasSeqData(chromosomes, inp_annotation, ref_genome, p1_dir, p2_dir, pas_samplenames, output_dir)

methods.Quantification(chromosomes, s1_dir, s2_dir, samplenames, filename, result_filename)
print("Total Program time: ",time.time() - startTime)

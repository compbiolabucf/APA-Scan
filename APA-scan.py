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

(speciesFlag, inputFlag, outFlag, pasFlag) = (0, 0, 0, 0)

startTime = time.time()

fasta_file_h = 'hg38.fa'
fasta_file_m = 'mm10.fa'
chromosomes_h = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21', 'chr22','chrX','chrY']
chromosomes_m = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']

if(len(sys.argv)<5):
	print("Please provide all of the mandatory arguments. Example: $ python3 apa_main.py -s h s1 s2")
	sys.exit()

if len(sys.argv)>=5:
	for ii in range(len(sys.argv)):
		if sys.argv[ii] == '-s' or sys.argv[ii] == '-S':
			speciesFlag = 1
			species = sys.argv[ii+1]
		elif sys.argv[ii] == '-o' or sys.argv[ii] == '-O':
			outFlag = 1
			output_dir = sys.argv[ii+1]+'/'
		elif sys.argv[ii] == '-p' or sys.argv[ii] == '-P':
			pasFlag = 1
			pasSeq_input1 = sys.argv[ii+1]
			pasSeq_input2 = sys.argv[ii+2]

input1_dir = sys.argv[3]
input2_dir = sys.argv[4]

if speciesFlag == 0:
	print("Please provide species name")
	sys.exit()
if outFlag == 0:
	output_dir = 'Output-APA-all/'


os.makedirs(output_dir, exist_ok=True)

s1_dir, s2_dir, samplenames, bamfile_names = methods.processInputFiles(input1_dir, input2_dir)

if len(bamfile_names) == 0:
	print("Bamfiles not found")
	sys.exit()


if species =='h':
	chromosomes = chromosomes_h
	inp_annotation = 'hg38_refseq_2018May1.txt'
	fasta_file = fasta_file_h
elif species == 'm':
	chromosomes = chromosomes_m
	inp_annotation = 'mm10_refseq_2019June20.txt'
	fasta_file = fasta_file_m


if pasFlag == 0:
	
	for i in range(2):
		if i==0:
			inputPath = s1_dir
		else:
			inputPath = s2_dir
	preprocess.SamtoText(input_path, samplenames[i], bamfile_names[i], chromosomes)
	
	filename = methods.Generate_withPasSeqSignal(chromosomes, inp_annotation, fasta_file, output_dir)
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
	
	filename = methods.Generate_withPasSeqData(chromosomes, inp_annotation, fasta_file, p1_dir, p2_dir, pas_samplenames, output_dir)

methods.Quantification(chromosomes, s1_dir, s2_dir, samplenames, filename, output_dir)
print("Total Program time: ",time.time() - startTime)
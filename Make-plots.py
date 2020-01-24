#!/usr/bin/env python3
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import bisect
from bisect import bisect_left
import csv
import time
import sys
import os
import methods

y_limit = 0
y_limit2 = 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def Generate_annotation_plot(ax, start, end, pos, strand):
	ax = ax or plt.gca()

	length = end-start+1
	p = (length*0.02)
	if strand == '+':
		rect1 = patches.Rectangle((start,0.25), length, 0.1, color = 'skyblue', fill = True)
		ax.add_patch(rect1)
		rect2 = patches.Rectangle((start,0.2), p, 0.2, color = 'skyblue', fill = True)
		ax.add_patch(rect2)

		length2 = pos-start+1
		rect1 = patches.Rectangle((start,0.5), length2, 0.1, color = 'skyblue', fill = True)
		ax.add_patch(rect1)
		rect2 = patches.Rectangle((start,0.45), p, 0.2, color = 'skyblue', fill = True)
		ax.add_patch(rect2)
	else:
		prect = start+(length*0.98)
		rect1 = patches.Rectangle((start,0.25), length, 0.1, color = 'skyblue', fill = True)
		ax.add_patch(rect1)
		rect2 = patches.Rectangle((prect,0.2), p, 0.2, color = 'skyblue', fill = True)
		ax.add_patch(rect2)

		length2 = end-pos+1
		rect1 = patches.Rectangle((pos,0.5), length2, 0.1, color = 'skyblue', fill = True)
		ax.add_patch(rect1)
		rect2 = patches.Rectangle((prect,0.45), p, 0.2, color = 'skyblue', fill = True)
		ax.add_patch(rect2)

	ax.set_xlim(start, end)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.set_yticklabels([])
	ax.set_xticklabels([])
	ax.tick_params(left=False, bottom=False, axis='both', which='major', labelsize=6)

	return

def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, start, end, pos, strand, number):
	bam_file_reader= open(pathin+sample+'/'+chrom+".txt", "rt")
	bam_read = csv.reader(bam_file_reader, delimiter="\t")
	bam_list = list(bam_read)
	position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]

	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)
	ax.tick_params(axis='both', which='major', labelsize=6)

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)
	if(int(bam_list[pos2][1]) != end):
		pos2 = pos2 - 1

	p = []
	c = []
	read = 0
	length = end - start + 1

	for t in range(length):
		p.append(t+start)
		c.append(0)
		
	for t in range(pos1, pos2+1):
		position = int(bam_list[t][1])
		read = int(bam_list[t][2])
		index = p.index(position)
		c[index] = read

	p = np.array(p)
	c = np.array(c)

	labelname = ""
	if number == 1:
		global y_limit
		m = max(c)
		if m > y_limit:
			y_limit = m

		labelname = sample+" (RNA-seq)"
		caption = ax.fill_between(p,c, color="skyblue", alpha=0.9, label = labelname)
		yy = y_limit

	elif number == 2:
		global y_limit2
		m = max(c)
		if m > y_limit2:
			y_limit2 = m

		labelname = sample+" (3'-end-seq)"
		caption = ax.fill_between(p,c, color="#884ea0", alpha=0.9, label = labelname)
		yy = y_limit2

	ax.legend(handles = [caption])
	ax.vlines(x=pos, ymin=0, ymax=m, colors='crimson', linestyles='solid', linewidth=1)
	ax.set_xlim(start, end)
	
	return yy

def Plot_Function_with_PAS(s1_dir, s2_dir, samplenames, p1_dir, p2_dir, pas_samplenames, file_list, chromDict, chrom, geneID, output_dir):

	geneList = chromDict[chrom]
	#print(GeneList)
	for (gene, strand) in geneList:
		if gene == geneID:
			df = pd.DataFrame(file_list)
			ann_tt = df.loc[(df[0] == chrom) & (df[1] == geneID)]
			for a_row in ann_tt.itertuples():
				start = int(a_row[4])
				end = int(a_row[5])
				pos = int(a_row[6])
				p_value = round(float(a_row[7]), 2)

				y_axis_height = 20
				title = ""+chrom+":"+geneID+":"+str(start)+"-"+str(end)+"("+strand+")"

				fig = plt.figure(figsize=(6.5,7))

				ax1 = fig.add_subplot(5,1,1)
				ax1.set_title(title, color = "black")
				ax1.set_ylabel('Counts')
				y_limit = Generate_read_coverate_plot(ax1, s1_dir, samplenames[0], chrom, geneID, int(start), int(end), int(pos), strand, 1)

				ax2 = fig.add_subplot(5,1,2)
				ax2.set_ylabel('Counts')
				y_limit2 = Generate_read_coverate_plot(ax2, p1_dir, pas_samplenames[0], chrom, geneID, int(start), int(end), int(pos), strand, 2)
				
				ax3 = fig.add_subplot(5,1,3)
				ax3.set_ylabel('Counts')
				y_limit = Generate_read_coverate_plot(ax3, s2_dir, samplenames[1], chrom, geneID, int(start), int(end), int(pos), strand, 1)

				ax4 = fig.add_subplot(5,1,4)
				ax4.set_xlabel('Position')
				ax4.set_ylabel('Counts')
				y_limit2 = Generate_read_coverate_plot(ax4, p2_dir, pas_samplenames[1], chrom, geneID, int(start), int(end), int(pos), strand, 2)

				# set ylim for both axes after getting max ylim values
				ax1.set_ylim(0, y_limit+y_limit*0.1)
				ax2.set_ylim(0, y_limit2+y_limit2*0.1)
				ax3.set_ylim(0, y_limit+y_limit*0.1)
				ax4.set_ylim(0, y_limit2+y_limit2*0.1)

				ax5 = fig.add_subplot(5,1,5)
				ax5.set_ylabel('Annotation')
				Generate_annotation_plot(ax5, start, end, pos, strand)

				plt.savefig(output_dir+title+'.png')
				plt.savefig(output_dir+title+'.eps', format = 'eps', dpi = 1000)
				y_limit = 0
				y_limit2 = 0



def Plot_Function_without_PAS(s1_dir, s2_dir, samplenames, file_list, chromDict, chrom, geneID, output_dir):

	geneList = chromDict[chrom]
	#print(GeneList)
	for (gene, strand) in geneList:
		if gene == geneID:
			df = pd.DataFrame(file_list)
			ann_tt = df.loc[(df[0] == chrom) & (df[1] == geneID)]
			for a_row in ann_tt.itertuples():
				start = int(a_row[4])
				end = int(a_row[5])
				pos = int(a_row[6])
				p_value = round(float(a_row[7]), 2)

				y_axis_height = 20
				title = ""+chrom+":"+str(start)+"-"+str(end)+"("+geneID+","+strand+")"

				fig = plt.figure()

				ax1 = fig.add_subplot(3,1,1)
				ax1.set_title(title, color = "black")
				ax1.set_ylabel('Counts')
				y_limit = Generate_read_coverate_plot(ax1, s1_dir, samplenames[0], chrom, geneID, int(start), int(end), int(pos), strand, 1)
				
				ax2 = fig.add_subplot(3,1,2)
				ax2.set_xlabel('Position')
				ax2.set_ylabel('Counts')
				y_limit = Generate_read_coverate_plot(ax2, s2_dir, samplenames[1], chrom, geneID, int(start), int(end), int(pos), strand, 1)

				# set ylim for both axes after getting max ylim values
				ax1.set_ylim(0, y_limit+y_limit*0.1)
				ax2.set_ylim(0, y_limit+y_limit*0.1)

				ax3 = fig.add_subplot(3,1,3)
				ax3.set_ylabel('Annotation')
				Generate_annotation_plot(ax3, start, end, pos, strand)

				plt.savefig(output_dir+title+'.png')
				plt.savefig(output_dir+title+'.eps', format = 'eps', dpi = 1000)
				y_limit = 0



######### Main starts here #################
startTime = time.time()
(speciesFlag, inputFlag, outFlag, pasFlag) = (0, 0, 0, 0)

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
	output_dir = 'Output-APA-Plots/'

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
    	if '_' not in rows[2] and rows[2]!='chrM':
	    	chromosomes.append(rows[2])
    chr_set = set(chromosomes)
    chromosomes = list(chr_set)

chromDict = {}
with open(inp_annotation, 'r') as f:
    reader = csv.reader(f, dialect='excel', delimiter='\t')
    headers = next(f)
    readerList = list(reader)
    for chrom in chromosomes:
    	geneList = []
    	df = pd.DataFrame(readerList)
    	rowsOfChr = df.loc[df[2] == chrom]
    	for row in rowsOfChr.itertuples():
    		geneID = row[13].strip()
    		strand = row[4]
    		if (geneID, strand) not in geneList:
    			geneList.append((geneID, strand))

    	chromDict[chrom] = geneList

region = input("Enter the range: Chrom:GeneID:RegionStart-RegionEnd ")
chrom, geneID, rng = region.split(':')
start, end = rng.split('-')

if pasFlag == 0:
	file_reader = open(output_dir+'Result.csv', "rt")
	headers = next(file_reader)
	file_read = csv.reader(file_reader, delimiter="\t")
	file_list = list(file_read)
	Plot_Function_without_PAS(s1_dir, s2_dir, samplenames, file_list, chromDict, chrom, geneID, output_dir)
else:
	p1_dir, p2_dir, pas_samplenames, pas_bamfile_names = methods.processInputFiles(pasSeq_input1, pasSeq_input2)
	print(output_dir)
	if(os.path.exists(output_dir+'Result_PAS.csv')==True):
		if len(pas_bamfile_names) == 0:
			print("Bamfiles not found")
			sys.exit()
		
		file_reader = open(output_dir+'Result_PAS.csv', "rt")
		headers = next(file_reader)
		file_read = csv.reader(file_reader, delimiter="\t")
		file_list = list(file_read)
		Plot_Function_with_PAS(s1_dir, s2_dir, samplenames, p1_dir, p2_dir, pas_samplenames, file_list, chromDict, chrom, geneID, output_dir)
	else:
		print("Result_PAS file not found")
		sys.exit()

totalTime = time.time() - startTime
print(geneID+" Plotted Successfully.")
print("Total program time is : ",totalTime)
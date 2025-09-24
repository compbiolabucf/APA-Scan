import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import bisect
from bisect import bisect_left
import csv
import methods
import time
import sys
import os, glob
import configparser
import xlsxwriter
from matplotlib import rcParams
plt.rc('legend',**{'fontsize':10})

y_limit, y_limit2 = 0, 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

def Generate_read_coverate_plot(ax, pathin, sample, labelname, chrom, geneID, start, end, pos, n1, p1):
	#print(os.path.join(pathin, sample, chrom+".txt"))
	bam_df = pd.read_csv(os.path.join(pathin, sample, chrom+".txt"), delimiter='\t')
	position_row = bam_df.iloc[:, 0].tolist()
	bam_list = bam_df.values.tolist()

	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)
	if pos2 >= len(bam_list) or int(bam_list[pos2][0]) != end:
		pos2 = pos2 - 1

	p = []
	c = []
	read = 0
	length = end - start + 1
	
	for t in range(length):
		p.append(t+start)
		c.append(0)
		
	for t in range(pos1, pos2+1):
		position = int(bam_list[t][0])
		read = int(bam_list[t][1])
		index = p.index(position)
		c[index] = read

	p = np.array(p)
	c = np.array(c)

	if p1 == 0:
		labelname = labelname+" (RNA-seq)"
		global y_limit
		m = max(c)
		if m > y_limit:
			y_limit = m
		yy = y_limit
	else:
		labelname = labelname+" (3'-end-seq)"
		global y_limit2
		m = max(c)
		if m > y_limit2:
			y_limit2 = m
		yy = y_limit2

	if n1 == 1:
		caption = ax.fill_between(p,c, color="midnightblue", alpha=0.9, label = labelname)
	elif n1 == 2:
		caption = ax.fill_between(p,c, color="#884ea0", alpha=0.9, label = labelname)

	ax.legend(handles = [caption])
	ax.vlines(x=pos, ymin=0, ymax=m, colors='crimson', linestyles='solid', linewidth=1)
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.set_xlim(start, end)
	ax.autoscale(enable = True)
	ax.set_xticklabels([])
	ax.tick_params(axis='both', bottom=False, which='major', labelsize=6)

	return y_limit, y_limit2


def Generate_annotation_plot(ax, ann_df, geneID, start, end, pos):
	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)

	strand = ann_df[ann_df['name2']==geneID]['strand'].iloc[0]

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
	ax.autoscale(enable = True)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.set_yticklabels([])
	ax.tick_params(left=False, axis='both', which='major', labelsize=6)

	return

def Plot_Function(pas_flag, region, input1_dir, input2_dir, s1_namelist, s2_namelist, pasSeq1_dir, pasSeq2_dir, p1_namelist, p2_namelist, ann_df, output_dir):
	g1_name = input1_dir.split("/")[-1]
	g2_name = input2_dir.split("/")[-1]

	chrom, geneID, rng, pos = region.split(':')
	start, end = rng.split('-')
	start, end, pos = int(start), int(end), int(pos)
	print("chrom gene start end position:", chrom, geneID, start, end, pos)
	title = geneID+":"+str(start)+"-"+str(end)+":"+str(pos)

	#fig = plt.figure(figsize=(8,8))
	x_inches = 6.4     # [mm]*constant
	y_inches = 4.8/5*(len(s1_namelist)*2+1)
	dpi = 100
	fig = plt.figure(1, figsize = (x_inches,y_inches), dpi = dpi, constrained_layout = True)
	ax = fig.add_subplot(111)
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
	ax.set_title(title, color = "black", fontsize = 16)
	ax.set_ylabel('Read Coverage', fontsize = 14)

	number_of_subplots = len(s1_namelist)+len(s2_namelist)+1
	print("Generating read coverage plots...")
	if pas_flag == 0:
		fig, axes = plt.subplots(nrows=number_of_subplots, ncols=1)
		for i in range (0,len(s1_namelist)):
			y_limit, y_limit2 = Generate_read_coverate_plot(axes[i], input1_dir, s1_namelist[i], g1_name, chrom, geneID, start, end, pos, 1, 0)
		for i in range(len(s1_namelist),number_of_subplots-1):
			j = i - len(s1_namelist)
			y_limit, y_limit2 = Generate_read_coverate_plot(axes[i], input2_dir, s2_namelist[j], g2_name, chrom, geneID, start, end, pos, 2, 0)
	elif pas_flag == 1:
		number_of_subplots = (len(s1_namelist)+len(s2_namelist))*2+1
		fig, axes = plt.subplots(nrows=number_of_subplots, ncols=1)
		for i in range(0,len(s1_namelist)):
			y_limit, y_limit2 = Generate_read_coverate_plot(axes[2*i], input1_dir, s1_namelist[i], g1_name, chrom, geneID, start, end, pos, 1, 0)
			y_limit, y_limit2 = Generate_read_coverate_plot(axes[2*i+1], pasSeq1_dir, p1_namelist[i], g1_name, chrom, geneID, start, end, pos, 1, 1)
			
		for i in range(len(s1_namelist),len(s1_namelist)+len(s2_namelist)):
			j = i - len(s1_namelist)
			y_limit, y_limit2 = Generate_read_coverate_plot(axes[2*i], input2_dir, s2_namelist[j], g2_name, chrom, geneID, start, end, pos, 2, 0)
			y_limit, y_limit2 = Generate_read_coverate_plot(axes[2*i+1], pasSeq2_dir, p2_namelist[j], g2_name, chrom, geneID, start, end, pos, 2, 1)
	
	print("Generating annotation plots...")
	ax3 = axes[number_of_subplots-1]
	Generate_annotation_plot(ax3, ann_df, geneID, start, end, pos)
	
	ax3.set_xlabel('Position', fontsize="14")
	ax3.set_ylabel('Annotation', fontsize="14")
	"""
	ax3.spines['top'].set_color('none')
	ax3.spines['bottom'].set_color('none')
	ax3.spines['left'].set_color('none')
	ax3.spines['right'].set_color('none')
	ax3.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
	"""

	for i in range(number_of_subplots-1):
		if pas_flag == 0:
			axes[i].set_ylim(0, y_limit*1.1)
		else:
			if i%2==0:
				axes[i].set_ylim(0, y_limit*1.1)
			elif i%2==1:
				axes[i].set_ylim(0, y_limit2*1.1)

	y_limit, y_limit2 = 0, 0
	os.makedirs(output_dir, exist_ok=True)
	plt.savefig(output_dir+title+'_specified.png')
	plt.savefig(output_dir+title+'_specified.eps', format = 'eps', dpi = 1000)
	print("Plotted successfully.")


######### Main starts here #################
startTime = time.time()

config = configparser.ConfigParser()
config.read('configuration.ini')

input1_dir = config['INPUT_RNAseq']['input1']
input2_dir = config['INPUT_RNAseq']['input2']
if input1_dir[-1] == "/":
	input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
	input2_dir = input2_dir[:-1]
pasSeq1_dir = config['INPUT_PASseq']['pas1']
pasSeq2_dir = config['INPUT_PASseq']['pas2']
if pasSeq1_dir[:-1] == "/":
	pasSeq1_dir = pasSeq1_dir[:-1]
if pasSeq2_dir[:-1] == "/":
	pasSeq2_dir = pasSeq2_dir[:-1]
output_dir = config['OUTPUT_FOLDER']['output_dir']
if output_dir[-1] != "/":
	output_dir += "/"

os.makedirs(output_dir, exist_ok=True)
inp_annotation = config['ANNOTATION']['annotation']
ref_genome = config['ANNOTATION']['genome']
print("RNA-seq input 1 dir:", input1_dir)
print("RNA-seq input 2 dir:", input2_dir)
print("3'-end-seq input 1 dir:", pasSeq1_dir)
print("3'-end-seq input 2 dir:", pasSeq2_dir)
print("Output Dir:", output_dir) 
print("Annotation:", inp_annotation, ref_genome, "\n\n")

pas_flag = 1
if pasSeq1_dir=='NULL' or pasSeq2_dir == 'NULL':
	pas_flag = 0

s1_namelist = list_dirs(input1_dir)
s2_namelist = list_dirs(input2_dir)
p1_namelist, p2_namelist = "", ""
if pas_flag == 1:
	p1_namelist = list_dirs(pasSeq1_dir)
	p2_namelist = list_dirs(pasSeq2_dir)

print("Loading list of chromosomes from the annotation...")
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

ann_df = pd.read_csv(inp_annotation, delimiter='\t', index_col=0)

region = input("Enter the range: (chr:gene:start-end:position): ")
#region = "chr4:Rpl22:152332259-152334082:152333201"
Plot_Function(pas_flag, region, input1_dir, input2_dir, s1_namelist, s2_namelist, pasSeq1_dir, pasSeq2_dir, p1_namelist, p2_namelist, ann_df, output_dir)

totalTime = time.time() - startTime
print("Total program time is : ",totalTime)
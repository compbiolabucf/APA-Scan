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

y_limit = 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def Generate_annotation_plot(ax, start, end, pos, strand, y_limit):
	ax = ax or plt.gca()

	length = end-start+1
	p = (length*0.02)
	if strand == '+':
		rect1 = patches.Rectangle((start,0.25), length, 0.1, color = 'skyblue', fill = True)
		ax.add_patch(rect1)
		rect2 = patches.Rectangle((start,0.2), p, 0.2, color = 'skyblue', fill = True)
		ax.add_patch(rect2)

		length2 = pos-start+1
		#ax.hlines(y=0.027, xmin=start, xmax=end, linewidth=1, color='slateblue', linestyle = '--')
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


	#ax.autoscale(enable = True)
	ax.set_xlim(start, end)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.set_yticklabels([])
	ax.set_xticklabels([])
	ax.tick_params(left=False, bottom=False, axis='both', which='major', labelsize=6)

	return

def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, start, end, pos, strand):
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

	#print(start, end, pos1, pos2)

	p = []
	c = []
	read = 0
	length = end - start + 1
	#print(length)
	for t in range(length):
		p.append(t+start)
		c.append(0)
		
	for t in range(pos1, pos2+1):
		position = int(bam_list[t][1])
		read = int(bam_list[t][2])
		index = p.index(position)
		c[index] = read

	#print("p, c : ",p,c)
	#print("st en pos: ", start, end, length, pos)
	p = np.array(p)
	c = np.array(c)

	global y_limit
	m = max(c)
	if m > y_limit:
		y_limit = m
	
	caption = ax.fill_between(p,c, color="skyblue", alpha=0.9, label = sample)
	ax.legend(handles = [caption])
	ax.vlines(x=pos, ymin=0, ymax=m, colors='crimson', linestyles='solid', linewidth=1)
	ax.set_xlim(start, end)
	#ax.autoscale(enable = True)
	
	return y_limit


def Take_user_inputs(pathin, samplenames, file_list, chromDict):
	region = input("Enter the range: ")
	print("You entered " + str(region))
	chrom, geneID, rng = region.split(':')
	start, end = rng.split('-')

	geneList = chromDict[chrom]
	#print(GeneList)
	for (gene, strand) in geneList:
		if gene == geneID:
			#print("######################")
			df = pd.DataFrame(file_list)
			#print(df[0])
			ann_tt = df.loc[(df[0] == chrom) & (df[1] == geneID)]
			#print(ann_tt)
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
				y_limit = Generate_read_coverate_plot(ax1, pathin, samplenames[0], chrom, geneID, int(start), int(end), int(pos), strand)
				
				ax2 = fig.add_subplot(3,1,2)
				ax2.set_xlabel('Position')
				ax2.set_ylabel('Counts')
				y_limit = Generate_read_coverate_plot(ax2, pathin, samplenames[1], chrom, geneID, int(start), int(end), int(pos), strand)

				# set ylim for both axes after getting max ylim values
				ax1.set_ylim(0, y_limit+y_limit*0.1)
				ax2.set_ylim(0, y_limit+y_limit*0.1)

				ax3 = fig.add_subplot(3,1,3)
				ax3.set_ylabel('Annotation')
				Generate_annotation_plot(ax3, start, end, pos, strand, y_limit)

				plt.savefig(title+'.png')
				#plt.savefig(plotout+title+'.eps', format = 'eps', dpi = 1000)
				y_limit = 0



######### Main starts here #################
startTime = time.time()

chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']

input_file = 'mm10.fa'
inputPath = '/home/naima/input/mouse_M-_M+/RNA-seq_bam/'
samplenames = ['Plus_M', 'Minus_M']

file_reader = open('Results.csv', "rt")
headers = next(file_reader)
file_read = csv.reader(file_reader, delimiter="\t")
file_list = list(file_read)


chromDict = {}
with open('mm10_refseq_2019June20.txt', 'r') as f:
    reader = csv.reader(f, dialect='excel', delimiter='\t')
    headers = next(f)
    readerList = list(reader)
    for chrom in chromosomes:
    	geneList = []
    	df = pd.DataFrame(readerList)
    	rowsOfChr = df.loc[df[1] == chrom]
    	for row in rowsOfChr.itertuples():
    		geneId = row[12].strip()
    		strand = row[3]
    		if (geneId, strand) not in geneList:
    			geneList.append((geneId, strand))

    	chromDict[chrom] = geneList

Take_user_inputs(inputPath, samplenames, file_list, chromDict)
totalTime = time.time() - startTime
print("Total program time is : ",totalTime)

"""
chr4:Rpl22:152332259-152334082
chr14:Rpl15:18267822-18269316

chr8:Prdx2:84973999-84974811
chr3:Snapin:90488025-90489593
chr11:Ddx5:106780355-106782256
chr13:Pfkp:6579873-6581592
chr14:Ctsb:63142231-63145923





"""
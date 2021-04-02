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
import os
import configparser
import xlsxwriter
from matplotlib import rcParams
plt.rc('legend',**{'fontsize':14})

rcParams.update({
    'font.family':'arial',
    })
y_limit = 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, startAll, endAll, number):
	bam_file_reader= open(pathin+'/'+sample+'/'+chrom+".txt", "rt")
	bam_read = csv.reader(bam_file_reader, delimiter="\t")
	bam_list = list(bam_read)
	position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]

	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)

	pos1 = bi_contains(position_row, startAll)
	pos2 = bi_contains(position_row, endAll)
	if(int(bam_list[pos2][1]) != endAll):
		pos2 = pos2 - 1

	p = []
	c = []
	read = 0
	length = endAll - startAll + 1
	
	for t in range(length):
		p.append(t+startAll)
		c.append(0)
		
	for t in range(pos1, pos2+1):
		position = int(bam_list[t][1])
		read = int(bam_list[t][2])
		index = p.index(position)
		c[index] = read

	p = np.array(p)
	c = np.array(c)

	#pos3 = bi_contains(p,start)
	#pos4 = bi_contains(p,end)

	global y_limit
	m = max(c)
	if m > y_limit:
		y_limit = m

	labelname = ""
	if number == 1:
		labelname = "Mock"
		caption = ax.fill_between(p,c, color="midnightblue", alpha=0.9, label = labelname)
	else:
		labelname = "SARS_CoV-2"
		caption = ax.fill_between(p,c, color="crimson", alpha=0.9, label = labelname)

	ax.legend(handles = [caption])
	#ax.fill_between(p[pos3:pos4+1],c[pos3:pos4+1], color="orange", alpha=0.9)
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.set_xticklabels([])
	ax.tick_params(axis='both', bottom=False, which='major', labelsize=10)

	return y_limit

def Generate_annotation_plot(ax, strand, isoforms, exonCountList, exonStartList, exonEndList, startAll, endAll, pos):
	#print(startAll, endAll)
	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)
	print("isoforms are: ",isoforms)

	rnge = endAll-startAll+1

	ystart = 0
	height = 3
	for i in range(isoforms):
		if i>=15:
			print("15 isoforms of this Gene is plotted.")
			break;
		else:
			ax.hlines(y=(ystart+ystart+height)/2, xmin=startAll, xmax=endAll, linewidth=1, color='skyblue', linestyle = '--')
			
			ecount, cat = exonCountList[i]
			ecount = int(ecount)
			stList = exonStartList[i]
			enList = exonEndList[i]

			if (cat == 'short'):
				colr = 'crimson'
			else:
				colr = 'skyblue'

			for p in range(ecount):
				ex_s = int(stList[p])
				ex_e = int(enList[p])
				width = int(enList[p]) - int(stList[p]) + 1

				ww = 0.1*width
				if (strand == '+') and (p == ecount-1):
					width = ex_e - ex_s + 1
					rect1 = patches.Rectangle((ex_s, ystart+1), width, 1, color = colr, alpha=0.9, fill = True)
					ax.add_patch(rect1)
					rect2 = patches.Rectangle((ex_s, ystart), ww, height, color = colr, alpha=0.9, fill = True)
					ax.add_patch(rect2)
				elif (strand == '-') and (p==0):
					width = ex_e - ex_s + 1
					rect1 = patches.Rectangle((ex_s, ystart+1), width, 1, color = colr, alpha=0.9, fill = True)
					ax.add_patch(rect1)
					wstart = ex_s+(0.9*width)
					rect2 = patches.Rectangle((wstart, ystart), ww, height, color = colr, alpha=0.9, fill = True)
					ax.add_patch(rect2)
				else:
					rect = patches.Rectangle((ex_s,ystart), width, height, color = 'skyblue', fill = True)
					ax.add_patch(rect)

			ystart +=5

	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.set_yticklabels([])
	ax.tick_params(left=False, axis='both', which='major', labelsize=10)

	return

def Plot_Function(input_dir, s1_namelist, s2_namelist, gene, output_dir, ann_list):
	df = pd.DataFrame(ann_list)
	ann_tt = df.loc[df[10]==gene]
	#if len(ann_tt) == 2:
	if len(ann_tt)==2 or len(ann_tt)==3 or len(ann_tt)==4:
		print("########", gene)
		exonStartList = {}
		exonEndList = {}
		exonCountList = {}

		isoforms = 0
		position = 0
		mini = 500000000
		maxi = 0
		tx_start_list, tx_end_list = [], []
		for a_row in ann_tt.itertuples():
			chrom = a_row[2]
			strand = a_row[3]
			tx_start = int(a_row[4])
			tx_end = int(a_row[5])
			tx_start_list.append(tx_start)
			tx_end_list.append(tx_end)
			cds_start = int(a_row[6])
			cds_end = int(a_row[7])
			#print("tx_start, tx_end, cds_start, cds_end: ", tx_start, tx_end, cds_start, cds_end)

			exonCount = int(a_row[8])
			category = a_row[12]
			exonCountList[isoforms] = (exonCount, category)
			exonStartList[isoforms] = a_row[9].split(',')
			exonEndList[isoforms] = a_row[10].split(',')
			isoforms+=1

			if strand == '+':
				if cds_end < mini:
					mini = cds_end
			elif strand == '-':
				if cds_start > maxi:
					maxi = cds_start

		if strand == '+':
			pos = mini
		else:
			pos = maxi
		all_start = min(np.array(tx_start_list))
		all_end = max(np.array(tx_end_list))

		title = ""+geneID+""

		#fig = plt.figure(figsize=(8,8))
		fig = plt.figure(figsize=(8,4))
		ax = fig.add_subplot(111)
		ax.spines['top'].set_color('none')
		ax.spines['bottom'].set_color('none')
		ax.spines['left'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
		ax.set_title(title, color = "black", fontsize = 20)
		ax.set_ylabel('Read Coverage', fontsize = 18)
		ax.set_xlabel('Position', fontsize = 18)

		#ax1 = fig.add_subplot(7,1,1)
		#y_limit = Generate_read_coverate_plot(ax1, input_dir, s1_namelist[0], chrom, geneID, all_start, all_end, 1)
		#ax2 = fig.add_subplot(7,1,2)
		#y_limit = Generate_read_coverate_plot(ax2, input_dir, s1_namelist[1], chrom, geneID, all_start, all_end, 1)
		#ax3 = fig.add_subplot(7,1,3)
		#y_limit = Generate_read_coverate_plot(ax3, input_dir, s1_namelist[2], chrom, geneID, all_start, all_end, 1)
		#ax4 = fig.add_subplot(7,1,4)
		#y_limit = Generate_read_coverate_plot(ax4, input_dir, s2_namelist[0], chrom, geneID, all_start, all_end, 2)
		#ax5 = fig.add_subplot(7,1,5)
		#y_limit = Generate_read_coverate_plot(ax5, input_dir, s2_namelist[1], chrom, geneID, all_start, all_end, 2)
		#ax6 = fig.add_subplot(7,1,6)
		#y_limit = Generate_read_coverate_plot(ax6, input_dir, s2_namelist[2], chrom, geneID, all_start, all_end, 2)
		#ax7 = fig.add_subplot(7,1,7)
		#ax7.set_ylabel('Annotation', fontsize="16")
		#Generate_annotation_plot(ax7, strand, isoforms, exonCountList, exonStartList, exonEndList, all_start, all_end, pos)

		ax1 = fig.add_subplot(3,1,1)
		y_limit = Generate_read_coverate_plot(ax1, input_dir, s1_namelist[0], chrom, geneID, all_start, all_end, 1)
		ax2 = fig.add_subplot(3,1,2)
		y_limit = Generate_read_coverate_plot(ax2, input_dir, s2_namelist[0], chrom, geneID, all_start, all_end, 2)
		ax3 = fig.add_subplot(3,1,3)
		ax3.set_ylabel('Annotation', fontsize="18")
		Generate_annotation_plot(ax3, strand, isoforms, exonCountList, exonStartList, exonEndList, all_start, all_end, pos)

		ax1.set_ylim(0, y_limit*1.2)
		ax2.set_ylim(0, y_limit*1.2)

		os.makedirs(output_dir, exist_ok=True)
		plt.savefig(output_dir+title+'_test.png')
		#plt.savefig(output_dir+title+'.eps', format = 'eps', dpi = 1000)
		print("Plotted successfully.")
		y_limit = 0
		#sys.exit()



######### Main starts here #################
startTime = time.time()

inp_annotation = "hg38_refseq_2018May1.txt" 
ann_file_reader= open(inp_annotation, "rt")
ann_read = csv.reader(ann_file_reader, delimiter="\t")
ann_list = list(ann_read)

#s1_namelist = ['GSM4462348', 'GSM4462349', 'GSM4462350']
#s2_namelist = ['GSM4462351', 'GSM4462352', 'GSM4462343']
#input_dir = "/home/naima/input/Covid19_data_human/tophat_hg38/"
input_dir = "/home/naima/input/Covid19_data_human/Series_5_merged/"
output_dir = "/home/naima/codes/TCBB_Covid19/APA_Scan/Plots_merged_sample/"
s1_namelist = ['Series5_Mock_hg38']
s2_namelist = ['Series5_SARS-CoV-2_hg38']
os.makedirs(output_dir, exist_ok=True)

reader1 = open('APA_Scan_Series5_A549_Mock_Vs_Series5_A549_SARS-CoV-2.csv', "rt")
read1 = csv.reader(reader1, delimiter="\t")
data_list = list(read1)[1:144]

#for region in region_list:
for line in data_list:
	row = line[0].split(',')
	chrom = row[0]
	geneID = row[1]
	strand = row[2]
	start = row[3]
	end = row[4]
	pos = row[5]
	#print(geneID)
	if geneID == 'PMEPA1':
		Plot_Function(input_dir, s1_namelist, s2_namelist, geneID, output_dir, ann_list)


#gene_list = [('TMEM201'), ('CNPY2')]

totalTime = time.time() - startTime
print("Total program time is : ",totalTime)

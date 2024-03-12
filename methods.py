import csv
import time
import sys
import math
from operator import itemgetter
import pandas as pd
from Bio import SeqIO
import re
import bisect
from bisect import bisect_left
from scipy.stats import chisquare
import peakutils
import numpy as np
import os

class Stack:
	def __init__(self):
		self.items = []

	def size(self):
		return len(self.items)

	def isEmpty(self):
		return self.items == []

	def push(self, val):
		self.items.append(val)

	def top(self):
		if self.isEmpty():
			return None
		else:
			return self.items[self.size()-1]

	def pop(self):
		if self.isEmpty():
			return None
		else:
			return self.items.pop()

def bi_contains(lst, item):
    return bisect_left(lst, item)

def SamtoText(input_path, bamfile, chromosomes):
	output_dir = os.path.join(input_path, bamfile[:-4])
	os.makedirs(output_dir, exist_ok=True)
	cwd = os.getcwd()
	cmd1 = os.path.join(cwd, "samtools")+" index "+os.path.join(input_path, bamfile)		# make samtools index bamfile.bam.bai
	os.system(cmd1)

	print(bamfile,"...")
	for chrom in chromosomes:
		cmd2 = os.path.join(cwd, "samtools")+" view -b "+os.path.join(input_path, bamfile)+" "+chrom+" -o "+os.path.join(output_dir, chrom+".bam")
		cmd3 = os.path.join(cwd, "samtools")+" pileup "+os.path.join(output_dir, chrom+".bam")+" | cut -f 2,4 > "+os.path.join(output_dir, chrom+".txt")   ### Need to use pileup, not mpileup
		command = cmd2+";"+cmd3
		os.system(command)
	return

def getFinalTargetRegion(inputList):
	n = len(inputList)
	inputList.sort(key = itemgetter(1), reverse = True)
	st = Stack()
	st.push(inputList[0])

	for i in range(1,n):
		stacktop = st.top()
		if inputList[i][0] < stacktop[0] and inputList[i][1] < stacktop[0]:
			st.push(inputList[i])
		elif inputList[i][0] < stacktop[0] and inputList[i][1] > stacktop[0] and inputList[i][1] < stacktop[1]:
			st.pop()
			st.push((inputList[i][0], stacktop[1]))
		elif inputList[i][1] == stacktop[1] and inputList[i][0] < stacktop[0]:
			st.pop()
			st.push(inputList[i])

	newList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		newList.append(stacktop)
		st.pop()

	return newList


def findAllOccurance(p, s):
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)

def CountReadCoverage(chrom, start, end, bam_list, position_row):
	totalCount = 0

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)

	if pos2 >= len(bam_list) or int(bam_list[pos2][0]) != end:
		pos2 = pos2 - 1

	if(pos1 < len(bam_list) and pos2 < len(bam_list)):
		for t in range(pos1, pos2+1):
			read = int(bam_list[t][1])
			totalCount += read
		
	return totalCount

def read_bamfiles(input_dir, sample, chrom):
	ss = time.time()
	bam_df = pd.read_csv(os.path.join(input_dir, sample, chrom+".txt"), delimiter='\t')
	position_row = bam_df.iloc[:, 0].tolist()
	bam_list = bam_df.values.tolist()
	return bam_list, position_row

def Generate_coverage(chrom, start, end, pos, strand, bam_list, position_row):
	if strand == '+':
		length = pos-start
		targetLength = end-pos
		RC = CountReadCoverage(chrom, start, pos, bam_list, position_row)
		targetRC = CountReadCoverage(chrom, pos+1, end, bam_list, position_row)
	else:
		targetLength = pos-start
		length = end-pos
		RC = CountReadCoverage(chrom, pos+1, end, bam_list, position_row)
		targetRC = CountReadCoverage(chrom, start, pos, bam_list, position_row)
	
	n = targetRC/targetLength
	N = RC/length
	return n, N

def getFinalTargetRegion(inputList):
	n = len(inputList)
	inputList.sort(key = itemgetter(1), reverse = True)
	st = Stack()
	st.push(inputList[0])

	for i in range(1,n):
		stacktop = st.top()
		if inputList[i][0] < stacktop[0] and inputList[i][1] < stacktop[0]:
			st.push(inputList[i])
		elif inputList[i][0] < stacktop[0] and inputList[i][1] > stacktop[0] and inputList[i][1] < stacktop[1]:
			st.pop()
			st.push((inputList[i][0], stacktop[1]))
		elif inputList[i][1] == stacktop[1] and inputList[i][0] < stacktop[0]:
			st.pop()
			st.push(inputList[i])

	newList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		newList.append(stacktop)
		st.pop()

	return newList


def makeSplittedList(position_row, bam_list, start, end):						
	p = []
	r = []

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)

	if pos2 >= len(bam_list) or int(bam_list[pos2][1]) != end:
		pos2 = pos2 - 1

	if(pos1 < len(bam_list) and pos2 < len(bam_list)):
		for t in range(pos1, pos2+1):
			p.append(int(bam_list[t][0]))
			r.append(int(bam_list[t][1]))

	return (p,r)

def calculatePeaksWithVallys(peakAreaPos, peakAreaRead, indexes, allStartPoint, allEndPoint):
	peakRangeList = []
	peakStartPoint = allStartPoint
	prevPeakValue = 0
	indexLen = len(indexes)

	# because we will compare two peaks together to consider the valley in between, we wont take the last peak index in the for loop
	for k in range(indexLen-1):
		selectedRange = peakAreaRead[indexes[k]:indexes[k+1]+1]
		valley = min(selectedRange)
		valleyPos = selectedRange.index(valley)
		actualValleyPos = peakAreaPos[indexes[k]+valleyPos]

		if valley < 0.3*min(peakAreaRead[indexes[k]],peakAreaRead[indexes[k+1]]):
			peakEndPoint = actualValleyPos
			peakRangeList.append((peakAreaPos[peakAreaRead.index(max(prevPeakValue,peakAreaRead[indexes[k]]))], peakStartPoint, peakEndPoint))
			prevPeakValue = 0
			peakStartPoint = actualValleyPos + 1
		else:
			prevPeakValue = max(prevPeakValue, peakAreaRead[indexes[k]])


	ind = [i for i, value in enumerate(peakAreaRead) if value == max(prevPeakValue, peakAreaRead[indexes[indexLen-1]]) and peakAreaPos[i] >= peakStartPoint]

	# because ind is a list with indexes of the maxpeak value, for that we will take ind[0] (ind[0] or ind[1] are the indices of the same value)
	peakRangeList.append((peakAreaPos[ind[0]], peakStartPoint, allEndPoint))
	return peakRangeList


def findPeakPosition(p, r):
	listOfPeakPositions = []
	listOfPeakRange = []
	flag = 0
	for i in range(len(p)):
		if int(r[i])>1:
			if flag == 0:
				peakStartPoint = p[i]
				peakAreaPos = []
				peakAreaRead = []
				flag = 1
			if flag == 1:
				peakAreaPos.append(p[i])
				peakAreaRead.append(r[i])

		if (int(r[i])==1 and flag == 1) or (i == len(p)-1 and flag == 1):
			if i == len(p)-1:
				peakEndPoint = p[i]
			else:
				peakEndPoint = p[i-1]

			np_peakAreaPos = np.array(peakAreaPos)
			np_peakAreaRead = np.array(peakAreaRead)

			indexes = peakutils.indexes(np_peakAreaRead, thres=6, min_dist=35, thres_abs = True)

			numberOfIndex = len(indexes)

			if numberOfIndex == 1:
				listOfPeakRange.append((peakAreaPos[indexes[0]], peakStartPoint, peakEndPoint))

			elif numberOfIndex > 1:
				peakRangeList = calculatePeaksWithVallys(peakAreaPos, peakAreaRead, indexes, peakStartPoint, peakEndPoint)
			
				listOfPeakRange.extend(peakRangeList)

			flag = 0

	return listOfPeakRange

def mergePeaksFromBothSamples(listOfPeakRange1, listOfPeakRange2, strand):
	len1 = len(listOfPeakRange1)
	len2 = len(listOfPeakRange2)

	cleavageSites = []
	peakFlag = []
	minVal = 10000

	if len1 == 0:
		for (peak, st, en) in listOfPeakRange2:
			if strand == '+':
				cleavageSites.append(en)
			else:
				cleavageSites.append(st)
		return cleavageSites
	elif len2 == 0:
		for (peak, st, en) in listOfPeakRange1:
			if strand == '+':
				cleavageSites.append(en)
			else:
				cleavageSites.append(st)
		return cleavageSites

	if len1<len2:
		row = len1
		col = len2
		firstList = listOfPeakRange1
		secondList = listOfPeakRange2
	else:
		row = len2
		col = len1
		firstList = listOfPeakRange2
		secondList = listOfPeakRange1

	peakFlag1 = {}
	peakFlag2 = {}
	distanceToSortList = []

	for i in range(row):
		for j in range(col):
			if strand == '+':
				distance = abs(int(firstList[i][2])-int(secondList[j][2]))
				distanceToSortList.append((int(firstList[i][2]), int(secondList[j][2]), distance))
				peakFlag1[firstList[i][2]] = 0
				peakFlag2[secondList[j][2]] = 0
			else:
				distance = abs(int(firstList[i][1])-int(secondList[j][1]))
				distanceToSortList.append((int(firstList[i][1]), int(secondList[j][1]), distance))
				peakFlag1[firstList[i][1]] = 0
				peakFlag2[secondList[j][1]] = 0

	distanceToSortList.sort(key = itemgetter(2), reverse = False)
	
	for (pos1, pos2, dist) in distanceToSortList:
		if peakFlag1[pos1] == 0 and peakFlag2[pos2] == 0:
			averageOfTwoPeaks = math.ceil((pos1+pos2)/2)
			cleavageSites.append(averageOfTwoPeaks)
			peakFlag1[pos1] = 1
			peakFlag2[pos2] = 1
		elif peakFlag1[pos1] != 0 and peakFlag2[pos2] == 0:
			cleavageSites.append(pos2)
			peakFlag2[pos2] = 1
		elif peakFlag1[pos1] == 0 and peakFlag2[pos2] != 0:
			cleavageSites.append(pos1)
			peakFlag1[pos1] = 1

	return cleavageSites


###################### Two methods for processing 3-end-seq data ####################
#done
def Get_Peak_Positions(chromosomes, ann_df, p1_dir, p2_dir, p1_name, p2_name, output_dir, extended):
	print("Getting list of all 3'-end-seq peaks")
	output_columns = ['Chrom', 'Gene', 'Strand', 'Start', 'End',  'Positions']
	writer_list = []
	filename = os.path.join(output_dir, 'Peak_positions.csv')

	for chrom in chromosomes:
		tt = time.time()
		
		bam_list1, position_row1 = read_bamfiles(p1_dir, p1_name, chrom)
		bam_list2, position_row2 = read_bamfiles(p2_dir, p2_name, chrom)

		geneList = list(set(ann_df[ann_df['chrom']==chrom]['name2'].tolist()))

		for gene in geneList:
			gene_rows = ann_df.loc[(ann_df['chrom'] == chrom) & (ann_df['name2'] == gene)]
			targetList = []
			for index, row in gene_rows.iterrows():
				if str(row['name'].strip()).startswith('NM_') or str(row['name'].strip()).startswith('NR_'):
					cdsStart, cdsEnd = int(row['cdsStart']), int(row['cdsEnd'])
					txStart, txEnd = int(row['txStart']), int(row['txEnd'])
					exonCount, strand = int(row['exonCount']), row['strand']
					exonStartList = row['exonStarts'].split(',')[0:-1]
					exonEndList = row['exonEnds'].split(',')[0:-1]

					if cdsStart < cdsEnd:
						for i in range(exonCount):
							exonStart = int(exonStartList[i])
							exonEnd = int(exonEndList[i])
							if strand == '+':
								if(exonStart<=cdsEnd and cdsEnd<=exonEnd):
									if (exonStart, txEnd) not in targetList:
										targetList.append((exonStart, txEnd))
							elif strand == '-':
								if(exonStart<=cdsStart and cdsStart<=exonEnd):
									if (txStart, exonEnd) not in targetList:
										targetList.append((txStart, exonEnd))

			if len(targetList) > 0:
				revisedTargetList = getFinalTargetRegion(targetList)

				if extended.upper() == 'YES':
					if strand == '+':
						st, en = revisedTargetList[-1]
						del revisedTargetList[-1]
						revisedTargetList.append((st, en+10000))
					elif strand == '-':
						st, en = revisedTargetList[0]
						del revisedTargetList[0]
						revisedTargetList.append((st-10000, en))

				for (st, en) in revisedTargetList:
					(p1, r1) = makeSplittedList(position_row1, bam_list1, st, en)
					(p2, r2) = makeSplittedList(position_row2, bam_list2, st, en)
					listOfPeakRange1 = findPeakPosition(p1, r1)
					listOfPeakRange2 = findPeakPosition(p2, r2)

					if len(listOfPeakRange1) > 0 or len(listOfPeakRange2) > 0:
						cleavageSites = mergePeaksFromBothSamples(listOfPeakRange1, listOfPeakRange2, strand)
						writer_list.append((chrom, gene, strand, st, en, cleavageSites))

		print("Chrom", chrom, "done in", round((time.time() - tt)/60, 2), "minutes")

	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(filename, sep='\t')

#done
def with_PA_peaks(chromosomes, s1_dir, s2_dir, g1_name, g2_name, output_dir, result_filename):
	df_p = pd.read_csv(os.path.join(output_dir, "Peak_positions.csv"), delimiter="\t")

	writer_list = []
	output_columns = ['Chrom', 'Gene Name', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1_name, 'n2: '+g2_name, 'N1: '+g1_name, 'N2: '+g2_name]
	
	position_row = []
	for chrom in chromosomes:
		ss = time.time()
		bam_list1, position_row1 = read_bamfiles(s1_dir, g1_name, chrom)
		bam_list2, position_row2 = read_bamfiles(s2_dir, g2_name, chrom)
		
		target_tt = df_p.loc[df_p['Chrom']==chrom]
		for index, row in target_tt.iterrows():
			gene = row['Gene'].strip()
			strand = row['Strand']
			start = int(row['Start'])
			end = int(row['End'])
			positionsList = list(map(int, row['Positions'].strip('[\' ]').split(',')))
			
			signi_p = 1.1
			ratio_diff = 'Nan'
			signi_ratio_diff = 'Nan'
			abs_ratio_diff = 'Nan'
			n1_f, N1_f, n2_f, N2_f = 0, 0, 0, 0
			targetRC1_f, targetRC2_f, RC1_f, RC2_f = 0, 0, 0, 0
			length_f, targetLength_f, signi_pos = 0, 0, 0
			flag = 0
			for pos in positionsList:
				length = end - start + 1
				targetLength = length
				
				if (strand == '+' and (pos-start)<=(length*0.85)) or (strand == '-' and (pos-start)>=(length*0.15)):
					if strand == '+':
						length = pos-start
						targetLength = end - pos
						RC1 = CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
						RC2 = CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
						targetRC1 = CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
						targetRC2 = CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
					else:
						targetLength = pos-start
						length = end - pos
						RC1 = CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
						RC2 = CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
						targetRC1 = CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
						targetRC2 = CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
					
					#print(targetRC1, targetRC2, targetLength, RC1, RC2, length)
					
					n1 = targetRC1/targetLength
					n2 = targetRC2/targetLength
					N1 = RC1/length
					N2 = RC2/length

					#print(n1, n2, N2, N2)

					if N1!=0 and N2!=0:
						ratio_diff = (n1/N1) - (n2/N2)
					else:
						ratio_diff = 0

					N1 = N1 + n1
					N2 = N2 + n2

					if N1!= 0 and N2!=0:
						P0 = (n1+n2)/(N1+N2)
						n10 = N1 * P0
						n20 = N2 * P0
						exp = [n10, N1-n10, n20, N2-n20]
						if 0 not in exp:
							flag = 1
							res = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
							if res[1] < signi_p:
								signi_p = res[1]
								signi_ratio_diff = ratio_diff
								abs_ratio_diff = abs(signi_ratio_diff)
								n1_f = n1
								N1_f = N1
								n2_f = n2
								N2_f = N2
								RC1_f = RC1
								RC2_f = RC2
								targetRC1_f = targetRC1
								targetRC2_f = targetRC2
								signi_pos = pos
								length_f = length
								targetLength_f = targetLength

			if flag == 1:			
				writer_list.append((chrom, gene, strand, start, end, signi_pos, signi_p, signi_ratio_diff, abs_ratio_diff, n1_f, n2_f, N1_f-n1_f, N2_f-n2_f))
			
		print("Chrom ", chrom, "done in ", round((time.time() - ss)/60, 2), "minutes")
		#sys.exit()

	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(os.path.join(output_dir, result_filename+".csv"), sep='\t')

	print("APA-Scan quantification done.")
	return

#done 
def with_PA_peaks_all(chromosomes, s1_dir, s2_dir, g1_name, g2_name, output_dir, result_filename):
	df_p = pd.read_csv(os.path.join(output_dir, "Peak_positions.csv"), delimiter="\t")

	writer_list = []
	output_columns = ['Chrom', 'Gene Name', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1_name, 'n2: '+g2_name, 'N1: '+g1_name, 'N2: '+g2_name]
	
	position_row = []
	for chrom in chromosomes:
		print(chrom)
		ss = time.time()
		bam_list1, position_row1 = read_bamfiles(s1_dir, g1_name, chrom)
		bam_list2, position_row2 = read_bamfiles(s2_dir, g2_name, chrom)
		
		target_tt = df_p.loc[df_p['Chrom']==chrom]
		for index, row in target_tt.iterrows():
			gene = row['Gene'].strip()
			strand = row['Strand']
			start = int(row['Start'])
			end = int(row['End'])
			positionsList = list(map(int, row['Positions'].strip('[\' ]').split(',')))

			for pos in positionsList:
				pos = int(pos.strip())
				length = end - start + 1
				targetLength = length
				
				if (strand == '+' and (pos-start)<=(length*0.85)) or (strand == '-' and (pos-start)>=(length*0.15)):
					if strand == '+':
						length = pos-start
						targetLength = end - pos
						RC1 = CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
						RC2 = CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
						targetRC1 = CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
						targetRC2 = CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
					else:
						targetLength = pos-start
						length = end - pos
						RC1 = CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
						RC2 = CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
						targetRC1 = CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
						targetRC2 = CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
					
					n1 = targetRC1/targetLength
					n2 = targetRC2/targetLength
					N1 = RC1/length
					N2 = RC2/length

					if N1!=0 and N2!=0:
						ratio_diff = (n1/N1) - (n2/N2)
					else:
						ratio_diff = 0

					N1 = N1 + n1
					N2 = N2 + n2

					if N1!= 0 and N2!=0:
						P0 = (n1+n2)/(N1+N2)
						n10 = N1 * P0
						n20 = N2 * P0
						exp = [n10, N1-n10, n20, N2-n20]
						if 0 not in exp:
							flag = 1
							res = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
							writer_list.append((chrom, gene, strand, start, end, pos, res[1], ratio_diff, abs(ratio_diff), n1, n2, N1-n1, N2-n2))

		print("Chrom", chrom, "done in ", roudn((time.time() - ss)/60, 2), "minutes")
	
	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(result_filename, sep='\t')
	print("APA-Scan quantification done.")
	return
######################################################################################



###################### Methods for processing RNA-seq data ###########################
#done
def Get_Signal_Positions(chromosomes, ann_df, ref_genome, output_dir, extended):
    filename = os.path.join(output_dir, 'Signal_positions.csv')
    writer_list = []
    output_columns = ['Chrom', 'Gene', 'Strand', 'Start', 'End', 'Positions']
    fasta_sequences = SeqIO.parse(open(ref_genome),'fasta')

    for fasta in fasta_sequences:
    	chrom, sequence = fasta.id, str(fasta.seq)
    	if chrom in chromosomes:
    		tt = time.time()
    		geneList = list(set(ann_df[ann_df['chrom']==chrom]['name2'].tolist()))

    		for gene in geneList:
    			gene_rows = ann_df.loc[(ann_df['chrom'] == chrom) & (ann_df['name2'] == gene)]
    			targetList = []
    			for index, row in gene_rows.iterrows():
    				if str(row['name'].strip()).startswith('NM_') or str(row['name'].strip()).startswith('NR_'):
    					cdsStart, cdsEnd = int(row['cdsStart']), int(row['cdsEnd'])
    					txStart, txEnd = int(row['txStart']), int(row['txEnd'])
    					exonCount, strand = int(row['exonCount']), row['strand']
    					exonStartList = row['exonStarts'].split(',')[0:-1]
    					exonEndList = row['exonEnds'].split(',')[0:-1]

    					if cdsStart < cdsEnd:
	    					if strand == '+':
	    						regStart = cdsEnd
	    						regEnd = txEnd
	    						for i in range(exonCount):
	    							exonStart = int(exonStartList[i])
	    							exonEnd = int(exonEndList[i])
	    							if(exonStart<=regStart and regStart<=exonEnd):
	    								if (exonStart, regEnd) not in targetList:
	    									targetList.append((exonStart, regEnd))
	    									break

	    					elif strand == '-':
	    						regStart = txStart
	    						regEnd = cdsStart
	    						for i in range(exonCount):
	    							exonStart = int(exonStartList[i])
	    							exonEnd = int(exonEndList[i])
	    							if(exonStart<=regEnd and regEnd<=exonEnd):
	    								if (regStart, exonEnd) not in targetList:
	    									targetList.append((regStart, exonEnd))
	    									break

	    		if len(targetList) > 0:				
    				revisedTargetList = getFinalTargetRegion(targetList)

    				if extended.upper() == 'YES':
    					if strand == '+':
    						st, en = revisedTargetList[-1]
    						del revisedTargetList[-1]
    						revisedTargetList.append((st, en+10000))
    					elif strand == '-':
    						st, en = revisedTargetList[0]
    						del revisedTargetList[0]
    						revisedTargetList.append((st-10000, en))

    				if strand == '+':
    					for (st, en) in revisedTargetList:
    						length = en - st + 1
    						seq = sequence[st-2: en-2].upper()
    						listPlus = []
    						for i in findAllOccurance('AATAAA', seq):
    							pos = i+5
    							if pos<=(length*0.85):
    								listPlus.append(st+pos)
    						for i in findAllOccurance('ATTAAA', seq):
    							pos = i+5
    							if pos<=(length*0.85):
    								listPlus.append(st+pos)

    						if len(listPlus)>0:
    							writer_list.append((chrom, gene, strand, st, en, list(set(listPlus))))

    				elif strand == '-':
    					for (st, en) in revisedTargetList:
    						length = en - st + 1
    						seq = sequence[st-2: en-2].upper()
    						listMinus = []
    						for pos in findAllOccurance('TTTATT', seq):
    							if pos>=(length*0.15):
    								listMinus.append(st+pos)
    						for pos in findAllOccurance('TTTAAT', seq):
    							if pos>=(length*0.15):
    								listMinus.append(st+pos)

    						if len(listMinus)>0:
    							writer_list.append((chrom, gene, strand, st, en, list(set(listMinus))))

    		print("Chrom",chrom, " done in ", round((time.time() - tt)/60, 2), "minutes")

    df_output = pd.DataFrame(writer_list, columns=output_columns)
    df_output.to_csv(filename, sep="\t")

#done
def with_PAS_signal(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1, g2, output_dir, result_filename):
	len1,len2 = len(s1_namelist), len(s2_namelist)

	df_p = pd.read_csv(os.path.join(output_dir, "Signal_positions.csv"), delimiter="\t")
	writer_list = []
	output_columns = ['Chrom', 'Gene', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1, 'n2: '+g2, 'N1: '+g1, 'N2: '+g2]

	position_row = []
	for chrom in chromosomes:
		ss = time.time()
		s1_bam_list, s2_bam_list, s1_position_row, s2_position_row = {}, {}, {}, {}
		for sample1 in s1_namelist:
			bam_list, position_row = read_bamfiles(input1_dir, sample1, chrom)
			s1_bam_list[sample1] = bam_list
			s1_position_row[sample1] = position_row
		for sample2 in s2_namelist:
			bam_list, position_row = read_bamfiles(input2_dir, sample2, chrom)
			s2_bam_list[sample2] = bam_list
			s2_position_row[sample2] = position_row

		selected_rows = df_p.loc[df_p['Chrom']==chrom]
		for index, row in selected_rows.iterrows():
			gene = row['Gene'].strip()
			strand = row['Strand']
			start = int(row['Start'])
			end = int(row['End'])
			positionsList = list(map(int, row['Positions'].strip('[\' ]').split(',')))
			#print("positionsList", positionsList)
			signi_p = 1.1
			ratio_diff = 'Nan'
			signi_ratio_diff = 'Nan'
			abs_ratio_diff = 'Nan'
			n1_f, N1_f, n2_f, N2_f, signi_pos = 0, 0, 0, 0, 0
			flag = 0

			for pos in positionsList:
				length = end - start + 1
				targetLength = length
				
				s1_n, s1_N, s2_n, s2_N = 0, 0, 0, 0
				for sample1 in s1_namelist:
					n, N = Generate_coverage(chrom, start, end, pos, strand, s1_bam_list[sample1], s1_position_row[sample1])
					s1_n += n
					s1_N += N
				for sample2 in s2_namelist:
					n, N = Generate_coverage(chrom, start, end, pos, strand, s2_bam_list[sample2], s2_position_row[sample2])
					s2_n += n
					s2_N += N

				n1 = s1_n/len1
				N1 = s1_N/len1
				n2 = s2_n/len2
				N2 = s2_N/len2
				#print("n1, N1, n2, N2", n1, N1, n2, N2)
				
				ratio_diff = 0

				N1 = N1 + n1
				N2 = N2 + n2

				if N1>0 and N2>0:
					ratio_diff = (n1/N1) - (n2/N2)
					P0 = (n1+n2)/(N1+N2)
					n10 = N1 * P0
					n20 = N2 * P0
					exp = [n10, N1-n10, n20, N2-n20]
					if 0 not in exp:
						flag = 1
						res, p_value = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
						#print("res, p_value", res, p_value)
						if p_value < signi_p:
							signi_p = p_value
							signi_ratio_diff = ratio_diff
							abs_ratio_diff = abs(signi_ratio_diff)
							n1_f = n1
							N1_f = N1
							n2_f = n2
							N2_f = N2
							signi_pos = pos
				
			if flag == 1:
				writer_list.append((chrom, gene, strand, start, end, signi_pos, signi_p, signi_ratio_diff, abs_ratio_diff, n1_f, n2_f, N1_f-n1_f, N2_f-n2_f))
		print("Chrom ", chrom, " done in ", round((time.time() - ss)/60, 2), "minutes")

	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(os.path.join(output_dir, result_filename+".csv"), sep='\t')
	print("APA-Scan quantification done.")
	return

#done
def with_PAS_signal_all(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1, g2, output_dir, result_filename):
	len1, len2 = len(s1_namelist), len(s2_namelist)

	df_p = pd.read_csv(os.path.join(output_dir, "Signal_positions.csv"), delimiter="\t")
	writer_list = []
	output_columns = ['Chrom', 'Gene', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1, 'n2: '+g2, 'N1: '+g1, 'N2: '+g2]

	position_row = []
	for chrom in chromosomes:
		ss = time.time()
		s1_bam_list, s2_bam_list, s1_position_row, s2_position_row = {}, {}, {}, {}
		for sample1 in s1_namelist:
			bam_list, position_row = read_bamfiles(input1_dir, sample1, chrom)
			s1_bam_list[sample1] = bam_list
			s1_position_row[sample1] = position_row
		for sample2 in s2_namelist:
			bam_list, position_row = read_bamfiles(input2_dir, sample2, chrom)
			s2_bam_list[sample2] = bam_list
			s2_position_row[sample2] = position_row

		selected_rows = df_p.loc[df_p['Chrom']==chrom]
		for index, row in selected_rows.iterrows():
			gene = row['Gene'].strip()
			strand = row['Strand']
			start = int(row['Start'])
			end = int(row['End'])
			positionsList = list(map(int, row['Positions'].strip('[\' ]').split(',')))

			for pos in positionsList:
				pos = int(pos.strip())
				length = end - start + 1
				targetLength = length
				
				s1_n, s1_N, s2_n, s2_N = 0, 0, 0, 0
				for sample1 in s1_namelist:
					n, N = Generate_coverage(chrom, start, end, pos, strand, s1_bam_list[sample1], s1_position_row[sample1])
					s1_n += n
					s1_N += N
				for sample2 in s2_namelist:
					n, N = Generate_coverage(chrom, start, end, pos, strand, s2_bam_list[sample2], s2_position_row[sample2])
					s2_n += n
					s2_N += N

				n1 = s1_n/len1
				N1 = s1_N/len1
				n2 = s2_n/len2
				N2 = s2_N/len2
				
				ratio_diff = 0

				N1 = N1 + n1
				N2 = N2 + n2

				if N1>0 and N2>0:
					ratio_diff = (n1/N1) - (n2/N2)
					P0 = (n1+n2)/(N1+N2)
					n10 = N1 * P0
					n20 = N2 * P0
					exp = [n10, N1-n10, n20, N2-n20]
					if 0 not in exp:
						flag = 1
						res, p_value = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
						writer_list.append((chrom, gene, strand, start, end, pos, p_value, ratio_diff, abs(ratio_diff), n1, n2, N1-n1, N2-n2))

		print("Chrom", chrom, "done in ",  round((time.time() - ss)/60, 2), "minutes")

	df_output = pd.DataFrame(writer_list, columns=output_columns)
	df_output.to_csv(os.path.join(output_dir, result_filename+".csv"), sep='\t')
	print("APA-Scan quantification done.")
	return

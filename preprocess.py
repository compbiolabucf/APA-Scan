import os
import time


def SamtoText(input_path, sample, bamfile, chromosomes):
	cmd1 = "./samtools index "+input_path+"/"+sample+"/"+bamfile		## make samtools index accepted_hits.bam.bai
	os.system(cmd1)
	print("Samtools index run completed.")

	for chrom in chromosomes:
		print("Start of ", chrom)
		tt = time.time()
		cmd2 = "./samtools view -b "+input_path+"/"+sample+"/"+bamfile+" "+chrom+" -o "+input_path+"/"+sample+'/'+chrom+".bam"
		cmd3 = "./samtools pileup "+input_path+"/"+sample+'/'+chrom+".bam | cut -f 1,2,4 > "+input_path+"/"+sample+'/'+chrom+".txt"    ### Need to use pileup, not mpileup
		command = cmd2+";"+cmd3
		os.system(command)
		print("Samtools Time: ", time.time()-tt)
	print("Input text file generation completed")
	return

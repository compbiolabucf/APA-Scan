import os
import time
import concurrent.futures

def process_chromosome(input_path, bamfile, chrom, cwd, output_dir):
    cmd2 = cwd + "/samtools view -b " + input_path + "/" + bamfile + " " + chrom + " -o " + output_dir + '/' + chrom + ".bam"
    cmd3 = cwd + "/samtools pileup " + output_dir + '/' + chrom + ".bam | cut -f 1,2,4 > " + output_dir + '/' + chrom + ".txt"
    command = cmd2 + ";" + cmd3
    os.system(command)

def SamtoText(input_path, bamfile, chromosomes, number_of_cores):
    cwd = (os.getcwd()).split("/")[0:-1]
    cwd = '/'.join(cwd)
    output_dir = input_path + '/' + bamfile[:-4] + '/'
    os.makedirs(output_dir, exist_ok=True)
    cmd1 = cwd + "/samtools index " + input_path + "/" + bamfile   # make samtools index bamfile.bam.bai
    os.system(cmd1)
    print("Aligning bam files for", bamfile, "...")

    with concurrent.futures.ThreadPoolExecutor(max_workers=int(number_of_cores)) as executor:
        future_to_chrom = {executor.submit(process_chromosome, input_path, bamfile, chrom, cwd, output_dir): chrom for chrom in chromosomes}
        for future in concurrent.futures.as_completed(future_to_chrom):
            chrom = future_to_chrom[future]
            future.result()
    return

import sys, os, subprocess

def copy_number_estimator(rDNA_chromosome):
        #change files to linux format
        process_cat_whole_genome = subprocess.Popen(['cat', whole_genome_ref], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open ('whole_genome_fixed.fasta', 'w') as whole_genome_fixedfile:
                subprocess.Popen(['tr', '-d', '\'\\r\''], stdin=process_cat_whole_genome.stdout, stdout=whole_genome_fixedfile, stderr=subprocess.PIPE).wait()
	#index reference files
        subprocess.Popen(['bwa', 'index', 'whole_genome_fixed.fasta'], stderr=subprocess.PIPE).wait()
	#map the reads
	devnull = open(os.devnull, 'w')
	with open ('whole_genome_samfile.sam', 'w') as whole_genome_samfile:
		subprocess.Popen(['bwa', 'mem', 'whole_genome_fixed.fasta', fwd_reads_file, rvs_reads_file], stdout=whole_genome_samfile, stderr=devnull).wait()
	#sam to bam
	with open ('whole_genome_bamfile.bam', 'w') as whole_genome_bamfile:
		subprocess.Popen(['samtools', 'view', '-bS', 'whole_genome_samfile.sam'], stdout=whole_genome_bamfile, stderr=subprocess.PIPE).wait()
	#sort bam
	with open ('whole_genome_bamsort.bam', 'w') as whole_genome_bamsort:
		subprocess.Popen(['samtools', 'sort', 'whole_genome_bamfile.bam'], stdout=whole_genome_bamsort, stderr=subprocess.PIPE).wait()
	#compute depth profile
	with open ('whole_genome_depth.txt', 'w') as whole_genome_depthfile:
		subprocess.Popen(['samtools', 'depth', '-aa', 'whole_genome_bamsort.bam'], stdout=whole_genome_depthfile, stderr=subprocess.PIPE).wait()
	#retain rDNA chromosome
	rDNA_chromosome = '$1 == "' + rDNA_chromosome + '" { print }'
	with open ('chrome_depth.txt', 'w') as chrome_depthfile:
		subprocess.Popen(['awk', rDNA_chromosome, 'whole_genome_depth.txt'], stdout=chrome_depthfile, stderr=subprocess.PIPE).wait()
	#retain only rDNA region
	rDNA_region = 'FNR >= ' + rDNA_start_locus + ' && FNR <= ' + rDNA_end_locus
	with open ('rDNA_depth.txt', 'w') as rDNA_depthfile:
		subprocess.Popen(['awk', rDNA_region, 'chrome_depth.txt'], stdout=rDNA_depthfile, stderr=subprocess.PIPE)
	#retain only non rDNA region
	nrDNA_region = 'FNR < ' + rDNA_start_locus + ' || FNR > ' + rDNA_end_locus
	with open ('nrDNA_depth.txt', 'w') as nrDNA_depthfile:
		subprocess.Popen(['awk', nrDNA_region, 'chrome_depth.txt'], stdout=nrDNA_depthfile, stderr=subprocess.PIPE).wait()
	with open ('rDNA_depth_sorted.txt', 'w') as rDNA_sortedfile:
		subprocess.Popen(['sort', '-nk', '3', 'rDNA_depth.txt'], stdout=rDNA_sortedfile, stderr=subprocess.PIPE)
	with open ('nrDNA_depth_sorted.txt', 'w') as nrDNA_sortedfile:
		subprocess.Popen(['sort', '-nk', '3', 'nrDNA_depth.txt'], stdout=nrDNA_sortedfile, stderr=subprocess.PIPE).wait()
	median = '{arr[NR]=$3} END {if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}'
	rDNA_median = subprocess.Popen(['awk', median, 'rDNA_depth_sorted.txt'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	rDNA_median = rDNA_median.stdout.readlines()
	rDNA_median = int(rDNA_median[0])
	nrDNA_median = subprocess.Popen(['awk', median, 'nrDNA_depth_sorted.txt'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	nrDNA_median = nrDNA_median.stdout.readlines()
	nrDNA_median = int(nrDNA_median[0])
	copy_number_estimate = int(round(rDNA_median/nrDNA_median))
	print('copy_number_estimate = {0}'.format(copy_number_estimate))

#main function
def main():
	copy_number_estimator(rDNA_chromosome)

if __name__ == "__main__":
        whole_genome_ref = sys.argv[1]		#whole genome fasta file
	rDNA_start_locus = sys.argv[2]		#start locus of rDNA in whole genomw fasta file
	rDNA_end_locus = sys.argv[3]           #end locus of rDNA in whole genomw fasta file
	fwd_reads_file = sys.argv[4]                        #fastq file of forward reads
	rvs_reads_file = sys.argv[5]                        #fastq file of forward reads
	rDNA_chromosome = sys.argv[6]                        #name of chromosome containing rDNA
        main()

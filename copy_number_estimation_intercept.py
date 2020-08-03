import sys, os, subprocess

def parser(parameters):
	info = []
	with open(parameters, 'r') as file_handle:                  #open vcf file for reading
		for line in file_handle:
			if line.startswith('#'):
				continue
			read = False
			char_list = ''                  #initialising character list for reading
			for index in range(len(line)):                  #loop to parse each lin
				if line[index] == ':':
					read = True                     #parser can read
				elif line[index].isspace():             #if parser encounter space
					if len(char_list) > 0:
						info.append(char_list)                  #inserting value in info array
					char_list = ''                  #initialising character list for reading
					continue
				elif read:
					char_list += line[index]        #add character to the character list
	return info

def copy_number_estimator(info):
	whole_genome_ref = info[0]		#whole genome fasta file
	rDNA_start_locus = int(info[1])		#start locus of rDNA in whole genomw fasta file
	rDNA_end_locus = int(info[2])		#end locus of rDNA in whole genomw fasta file
	rDNA_buffer = int(info[3])		#buffer between rDNA and non rDNA
	rDNA_first_coord = int(info[4])		#first coordinate of rDNA in rDNA fasta file
	rDNA_last_coord = int(info[5])		#last coordinate of rDNA in rDNA fasta file
	rDNA_ref = info[6]			#fasta file of rDNA only
	fwd_reads_file = info[7]			#fastq file of forward reads
	rvs_reads_file = info[8]			#fastq file of reverse reads
	window_size = int(info[9])		#window size for computing average read depth
	step_size = int(info[10])		#step size for computing average read depth
        #change files to linux format
	process_cat_rDNA = subprocess.Popen(['cat', rDNA_ref], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	with open ('rDNA_fixed.fasta', 'w') as rDNA_fixedfile:
		subprocess.Popen(['tr', '-d', '\'\\r\''], stdin=process_cat_rDNA.stdout, stdout=rDNA_fixedfile, stderr=subprocess.PIPE)
	process_cat_whole_genome = subprocess.Popen(['cat', whole_genome_ref], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	with open ('whole_genome_fixed.fasta', 'w') as whole_genome_fixedfile:
                subprocess.Popen(['tr', '-d', '\'\\r\''], stdin=process_cat_whole_genome.stdout, stdout=whole_genome_fixedfile, stderr=subprocess.PIPE)
	#find maximum read length
	max_read_length = subprocess.Popen(['awk', '{if(NR%4==2) {} max_read_length = (max_read_length > length) ? max_read_length : length} END {print max_read_length}', fwd_reads_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	max_read_length = max_read_length.stdout.readlines()
	max_read_length_fwd = int(max_read_length[0])   #maximum read length in forward read file
	max_read_length = subprocess.Popen(['awk', '{if(NR%4==2) {} max_read_length = (max_read_length > length) ? max_read_length : length} END {print max_read_length}', rvs_reads_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	max_read_length = max_read_length.stdout.readlines()
	max_read_length_rvs = int(max_read_length[0])   #maximum read length in reverse read file
	max_read_length = max(max_read_length_fwd, max_read_length_rvs) #greater of forward and reverse reads
	shoulder_size = max_read_length - 1
	rDNA_first_coord = rDNA_first_coord - shoulder_size
	rDNA_last_coord = rDNA_last_coord + shoulder_size
	#trim the rDNA fasta file
	subprocess.Popen([python, 'fasta_trimmer.py', 'rDNA_fixed.fasta', 'rDNA_trm.fasta', '{0}'.format(rDNA_first_coord), '{0}'.format(rDNA_last_coord)]).wait()
	#index reference files
	subprocess.Popen(['bwa', 'index', 'rDNA_trm.fasta'], stderr=subprocess.PIPE)
	subprocess.Popen(['bwa', 'index', 'whole_genome_fixed.fasta'], stderr=subprocess.PIPE).wait()
	#map the reads
	devnull = open(os.devnull, 'w')
	with open ('rDNA_samfile.sam', 'w') as rDNA_samfile:
		subprocess.Popen(['bwa', 'mem', 'rDNA_trm.fasta', fwd_reads_file, rvs_reads_file], stdout=rDNA_samfile, stderr=devnull)	#to suppress screen output
	with open ('whole_genome_samfile.sam', 'w') as whole_genome_samfile:
		subprocess.Popen(['bwa', 'mem', 'whole_genome_fixed.fasta', fwd_reads_file, rvs_reads_file], stdout=whole_genome_samfile, stderr=devnull).wait()
	#sam to bam
	with open ('rDNA_bamfile.bam', 'w') as rDNA_bamfile:
		subprocess.Popen(['samtools', 'view', '-bS', 'rDNA_samfile.sam'], stdout=rDNA_bamfile, stderr=subprocess.PIPE)
	with open ('whole_genome_bamfile.bam', 'w') as whole_genome_bamfile:
		subprocess.Popen(['samtools', 'view', '-bS', 'whole_genome_samfile.sam'], stdout=whole_genome_bamfile, stderr=subprocess.PIPE).wait()
	#sort bam
	with open ('rDNA_bamsort.bam', 'w') as rDNA_bamsort:
		subprocess.Popen(['samtools', 'sort', 'rDNA_bamfile.bam'], stdout=rDNA_bamsort, stderr=subprocess.PIPE)
	with open ('whole_genome_bamsort.bam', 'w') as whole_genome_bamsort:
		subprocess.Popen(['samtools', 'sort', 'whole_genome_bamfile.bam'], stdout=whole_genome_bamsort, stderr=subprocess.PIPE).wait()
	#compute depth profile
	with open ('rDNA_depth.txt', 'w') as rDNA_depthfile:
		subprocess.Popen(['samtools', 'depth', '-aa', 'rDNA_bamsort.bam'], stdout=rDNA_depthfile, stderr=subprocess.PIPE)
	with open ('whole_genome_depth.txt', 'w') as whole_genome_depthfile:
		subprocess.Popen(['samtools', 'depth', '-aa', 'whole_genome_bamsort.bam'], stdout=whole_genome_depthfile, stderr=subprocess.PIPE).wait()
	#copy number estimate
	copy_number_estimate = subprocess.Popen([python, 'copyno.py', 'whole_genome_depth.txt', '{0}'.format(rDNA_start_locus), '{0}'.format(rDNA_end_locus), '{0}'.format(rDNA_buffer), 'rDNA_depth.txt', '{0}'.format(shoulder_size), '{0}'.format(window_size), '{0}'.format(step_size)], stdout=subprocess.PIPE)
	copy_number_estimate = copy_number_estimate.stdout.readlines()
	copy_number_estimate = int(copy_number_estimate[0])
	print(copy_number_estimate)

#main function
def main():
	info = parser(parameters)
	copy_number_estimator(info)


if __name__ == "__main__":
	python = 'python'+str(sys.version_info.major)   #python version
	parameters = sys.argv[1]                #input parameter file
	main()


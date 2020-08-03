import getopt, sys, os
import subprocess
from os import path

def variant_discoverer(first_coordinate, last_coordinate, max_read_len, rDNA_unit_size):
	#change file to linux format
	process_cat = subprocess.Popen(['cat', ref_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	with open ('fixed.fasta', 'w') as fixedfile:
		subprocess.Popen(['tr', '-d', '\'\\r\''], stdin=process_cat.stdout, stdout=fixedfile, stderr=subprocess.PIPE).wait()
	#trim reference file
	subprocess.Popen([python, 'fasta_trimmer.py', 'fixed.fasta', 'ref_trm.fasta', '{0}'.format(first_coordinate), '{0}'.format(last_coordinate)]).wait()
	#Index fasta with bwa
	try:
		subprocess.Popen(['bwa', 'index', 'ref_trm.fasta'], stderr=subprocess.PIPE).wait()
	except OSError:
		print('Error10: Cannot find bwa.')
		sys.exit(2)
	#Index fasta with samtool
	try:
		subprocess.Popen(['samtools', 'faidx', 'ref_trm.fasta'], stderr=subprocess.PIPE)
	except OSError:
		print('Error11: Cannot find samtools.')
		sys.exit(2)
	#Map the reads with bwa
	process_align_reads = subprocess.Popen(['bwa', 'mem', 'ref_trm.fasta',  fwd_reads_file,  rvs_reads_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#filter unmapped reads with samtools
	process_filter_unmapped_reads = subprocess.Popen(['samtools', 'view', '-h', '-F', '4'], stdin=process_align_reads.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#filter reads with soft clippings with software available at https://github.com/tseemann/samclip
	try:
		process_filter_softclipped_reads = subprocess.Popen(['samclip', '--ref', 'ref_trm.fasta', '--max', '0'], stdin=process_filter_unmapped_reads.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		print('Error12: Cannot find samclip.')
		sys.exit(2)
	#Convert file to BAM format for the next tool
	with open ('bamfile.bam', 'w') as bamfile:
		subprocess.Popen(['samtools', 'view', '-bS'], stdin=process_filter_softclipped_reads.stdout, stdout=bamfile, stderr=subprocess.PIPE).wait()
	#Use Prithika's new Python tool to convert CIGAR strings to CIGAR2 strings software available from Prithika
	subprocess.Popen([python, 'fat-cigar.py', 'linear', 'bamfile.bam', 'CG2.bam'], stderr=subprocess.PIPE).wait()
	#Convert the BAM files back to SAM format for CIGAR2 filtering
	with open ('CG2.sam', 'w') as samfile:
		subprocess.Popen(['samtools', 'view', '-h', 'CG2.bam'], stdout=samfile, stderr=subprocess.PIPE).wait()
	#Copy header from the sam file
	with open ('ch.sam', 'w') as chfile:
		subprocess.Popen(['samtools', 'view', '-H', 'CG2.sam'], stdout=chfile, stderr=subprocess.PIPE)
	#CIGAR2 filtering - only keep reads if they they have at least 20 nucleotide matches at each end of the read
	with open ('filtered.sam', 'w') as filteredfile:
		subprocess.Popen(['awk', '($6 ~ /^[2-9][0-9]=.*[2-9][0-9]=$/) || ($6 ~ /^[2-9][0-9]=.*1[0-9]{2}=$/) || ($6 ~ /^1[0-9]{2}=.*[2-9][0-9]=$/) || ($6 ~ /^[1-9][0-9]{1,2}=$/)', 'CG2.sam'], stdout=filteredfile, stderr=subprocess.PIPE).wait()
	#Add the header to the CIGAR2-filtered files as the filtering step lost it
	process_add_header = subprocess.Popen(['cat', 'ch.sam', 'filtered.sam'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#Convert the CIGAR2-filtered SAM files to BAM format
	with open ('back.bam', 'w') as bamfile:
		if float(sub_sample) < 1.0:
			subprocess.Popen(['samtools', 'view', '-bS', '-s', sub_sample], stdin=process_add_header.stdout, stdout=bamfile, stderr=subprocess.PIPE).wait()
		else:
			subprocess.Popen(['samtools', 'view', '-bS'], stdin=process_add_header.stdout, stdout=bamfile, stderr=subprocess.PIPE).wait()
	#Stack the filtered reads under the reference - up to maximum depth using freebayes
	with open ('vcfile.vcf', 'w') as vcfile:
		try:
			subprocess.Popen(['freebayes', '-f', 'ref_trm.fasta', '-F', '0.001', '--pooled-continuous', 'back.bam'], stdout=vcfile, stderr=subprocess.PIPE).wait()
		except OSError:
			print('Error13: Cannot find freebayes.')
			sys.exit(2)
	#Break multi-allelic sites
	try:
		process_breakmulti = subprocess.Popen(['vcfbreakmulti', 'vcfile.vcf'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		print('Error14: Cannot find vcflib.')
		sys.exit(2)
	#Decompose biallelic sites
	process_allelicprimitives = subprocess.Popen(['vcfallelicprimitives', '-kg'], stdin=process_breakmulti.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#Fix the field of alternate allelic depth
	process_fix_depth = subprocess.Popen(['./vcf_ad_fix.sh'], stdin=process_allelicprimitives.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#normalise the vcf
	try:
		process_normalise = subprocess.Popen(['bcftools', 'norm', '-f', 'ref_trm.fasta', '-m-'], stdin=process_fix_depth.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		print('Error15: Cannot find bcftools.')
		sys.exit(2)
	#relocate variants
	subprocess.Popen([python, 'variant_relocator.py', 'relocated.vcf', str(max_read_len), str(rDNA_unit_size+max_read_len-1), '0'], stdin=process_normalise.stdout, stderr=subprocess.PIPE).wait()
	#sort variants
	with open ('sorted.vcf', 'w') as sortedfile:                      #file object to write on the file
		subprocess.Popen(['bcftools', 'sort', 'relocated.vcf'], stdout=sortedfile, stderr=subprocess.PIPE).wait()
	#unify decomposed variants
	subprocess.Popen([python, 'variant_unifier.py', 'sorted.vcf', 'unified.vcf', '1']).wait()
	#unify relocated variants
	subprocess.Popen([python, 'variant_unifier.py', 'unified.vcf', 'unified2.vcf', '0']).wait()
	#filter variants
	with open ('filtered.vcf', 'w') as filteredfile:                #file object to write on the file
		subprocess.Popen(['vcffilter', '-f', 'AO > 1 & AO / DP > 0.005', 'unified2.vcf'], stdout=filteredfile, stderr=subprocess.PIPE) #filtering the variants
	#retrieve frquency
	subprocess.Popen([python, 'frequency_retriever.py', str(max_read_len-1), str(rDNA_unit_size), out_vcf]).wait()

#function to delete all intermediate files
def delete_all_interim_files():
	#delete interim reference file
	subprocess.Popen(['rm', 'fixed.fasta'])
	#deleting the fat-cigar bam file
	subprocess.Popen(['rm', 'CG2.bam'])
	#deleting the fat-cigar sam file
	subprocess.Popen(['rm', 'CG2.sam'])
	#delete header of sam
	subprocess.Popen(['rm', 'ch.sam'])
	#delete sam without header
	subprocess.Popen(['rm', 'filtered.sam'])
	#delete fat-cigar filtered bam file
	subprocess.Popen(['rm', 'back.bam'])
	#delete vcf file
	subprocess.Popen(['rm', 'vcfile.vcf'])
	#delete relocated vcf file
	subprocess.Popen(['rm', 'relocated.vcf'])
	#delete sorted vcf file
	subprocess.Popen(['rm', 'sorted.vcf'])
	#delete unified vcf file
	subprocess.Popen(['rm', 'unified.vcf'])
	#delete unified2 vcf
	subprocess.Popen(['rm', 'unified2.vcf'])
        #delete bam file
	subprocess.Popen(['rm', 'bamfile.bam'])
	#delete sorted bam file
	subprocess.Popen(['rm', 'bamfile_sorted.bam'])
	#delete index of sorted bam file
	subprocess.Popen(['rm', 'bamfile_sorted.bam.bai'])
	#delete head removed vcf
	subprocess.Popen(['rm', 'vcf_hr.vcf'])
	#delete the trimmed reference fasta file
	subprocess.Popen(['rm', 'ref_trm.fasta'])
	#delete the amb index of reference fasta file
	subprocess.Popen(['rm', 'ref_trm.fasta.amb'])
	#delete the ann index of reference fasta file
	subprocess.Popen(['rm', 'ref_trm.fasta.ann'])
	#delete the bwt index of reference fasta file
	subprocess.Popen(['rm', 'ref_trm.fasta.bwt'])
	#delete the fai index of reference fasta file
	subprocess.Popen(['rm', 'ref_trm.fasta.fai'])
	#delete the pac index of reference fasta file
	subprocess.Popen(['rm', 'ref_trm.fasta.pac'])
	#delete the sa index of reference fasta file
	subprocess.Popen(['rm', 'ref_trm.fasta.sa'])
	#delete relocated frequency retrieved vcf
	subprocess.Popen(['rm', 'relocated_fr.vcf'])
	#delete sorted frequency retrieved vcf
	subprocess.Popen(['rm', 'sorted_fr.vcf'])
	#delete unified frequency retrieved vcf
	subprocess.Popen(['rm', 'unified_fr.vcf'])
	#delete unified frequency retrieved vcf
	subprocess.Popen(['rm', 'unified_fr2.vcf'])
	#delete filtered vcf file
	subprocess.Popen(['rm', 'filtered.vcf'])
	#delete filtered frequency vcf file
	subprocess.Popen(['rm', 'filtered_fr.vcf'])

#function to find maximum read length
def find_max_read_length():
	max_read_length = subprocess.Popen(['awk', '{if(NR%4==2) {} max_read_length = (max_read_length > length) ? max_read_length : length} END {print max_read_length}', fwd_reads_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	max_read_length = max_read_length.stdout.readlines()
	max_read_length_fwd = int(max_read_length[0])   #maximum read length in forward read file
	max_read_length = subprocess.Popen(['awk', '{if(NR%4==2) {} max_read_length = (max_read_length > length) ? max_read_length : length} END {print max_read_length}', rvs_reads_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	max_read_length = max_read_length.stdout.readlines()
	max_read_length_rvs = int(max_read_length[0])   #maximum read length in reverse read file
	max_read_length = max(max_read_length_fwd, max_read_length_rvs) #greater of forward and reverse reads
	return max_read_length
					
#main function
def main(first_coordinate, last_coordinate):
	max_read_len = find_max_read_length()		#function to find maximum read length
	rDNA_unit_size = last_coordinate - first_coordinate + 1  #size of rDNA unit
	last_coordinate += max_read_len - 1	#last coordinate of trimmed fasta file
	first_coordinate -= max_read_len - 1	#first coordinate of trimmed fasta file
	variant_discoverer(first_coordinate, last_coordinate, max_read_len,  rDNA_unit_size)	#fuction to discover the variants
	delete_all_interim_files()	#function to delete all intermediate files

class take_command_line_args():
	def __init__(self, args, args_len):
		global ref_file, first_coordinate, last_coordinate, fwd_reads_file, rvs_reads_file, sub_sample, out_vcf
		ref = False
		first = False
		last = False
		fwd = False
		rvs = False
		out = False
		sub = False
		ref_file = ''
		fwd_reads_file = ''
		rvs_reads_file = ''
		sub_sample = '1.0'			#default value
		for i in range(0, args_len-1, 2):
			if args[i] == '--reference' or args[i] == '-f':
				ref_file = args[i+1]			#reference file
				ref = True
			elif args[i] == '--unit_start' or args[i] == '-u':
				if args[i+1].isdigit():
					first_coordinate = int(args[i+1])    #first coordinate of rDNA unit
					first = True
				else:
					print('Error16: First coordinate must be numeric')
					sys.exit(2)
			elif args[i] == '--unit_end' or args[i] == '-v':
				if args[i+1].isdigit():
					last_coordinate = int(args[i+1])     #last coordinate of rDNA unit
					last = True
				else:
					print('Error17: Last coordinate must be numeric')
					sys.exit(2)
			elif args[i] == '--reads1' or args[i] == '-r':
				fwd_reads_file = args[i+1]      #forward reads file
				fwd = True
			elif args[i] == '--reads2' or args[i] == '-q':
				rvs_reads_file = args[i+1]      #reverse reads file
				rvs = True
			elif args[i] == '--sub_sample' or args[i] == '-s':
				sub_sample = args[i+1]		#ratio of sampled reads in the bam file
				sub = True
			elif args[i] == '--out' or args[i] == '-o':
				out_vcf = args[i+1]             #output vcf file
				out = True
			else:
				print('Error0: option {0} is not recognised'.format(args[i]))
				sys.exit(2)
		#error ans warning list
		if not ref:
			print('Error1: Reference file name is missing')
		if not first:
			print('Error2: First coordinate of rDNA unit is missing')
		if not last:
			print('Error3: Last coordinate of rDNA unit is missing')
		if not fwd:
			print('Error4: Forward read file name is missing')
		if not rvs:
			print('Error5: Reverse read file name is missing')
		if not out:
			print('Error6: Output file name is missing')
		if ref and not path.exists(ref_file):
			print('Error7: Cannot find reference file')
		if fwd and not path.exists(fwd_reads_file):
			print('Error8: Cannot find read1 file')
		if rvs and not path.exists(rvs_reads_file):
			print('Error9: Cannot find read2 file')
		if not ref or not first or not last or not fwd or not rvs or not out or not path.exists(ref_file) or not path.exists(fwd_reads_file) or not path.exists(rvs_reads_file):
			sys.exit(2)
		if sub and float(sub_sample) >= 1.0 or float(sub_sample) <= 0.0:			#sampling ratio should always be less than 1
			print('warning: Sub_sample ratio must be greater than 0 and less than 1. By default reads will not be sampled.')
			sub_sample = 1.0
	def execute(self):
		main(first_coordinate, last_coordinate)		#adding two arguments because these are not acting as global variables for unknown reasons

if __name__ == "__main__":
	python = 'python'+str(sys.version_info.major)   #python version
	full_cmd_arguments = sys.argv	# Get full command-line arguments
	argument_list = full_cmd_arguments[1:]	# Keep all but the first (because first is the function name)
	short_options = "f:u:v:r:q:s:o:"
	long_options = ["reference=", "unit_start=", "unit_end=", "reads1=", "reads2=", "sub_sample=", "out="]
	try:
		arguments, values = getopt.getopt(argument_list, short_options, long_options)
	except getopt.error as err:
		# Output error, and return with an error code
		print (str(err))
		sys.exit(2)
	args_len = len(sys.argv)
	args = []
	for i in range(1, args_len):
		args.append(sys.argv[i])
	take_command_line_args(args, args_len).execute()

import sys, os, subprocess
from os import path

#function to parse configuration file
def parser():
	info = []
	with open(config, 'r') as file_handle:                  #open vcf file for reading
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
						if(char_list.isdigit()):			#if all characters in string are digits
							char_list = int(char_list)
						elif(char_list.replace('.','',1).isdigit()):	#if the string represents a floating point number
							char_list = float(char_list)
						info.append(char_list)                  #inserting value in info array
					char_list = ''                  #initialising character list for reading
					continue
				elif read:
                                        char_list += line[index]        #add character to the character list
	return info


def variant_generator(info):
	nr_of_repeats = info[0]		#number of rDNA repeats
	if not str(nr_of_repeats).isdigit():
		print('Error0: Number of rDNA repeats must be a number.')
		sys.exit(2)
	nr_of_noncoding_reg = info[1]	#number of non coding regions
	if not str(nr_of_noncoding_reg).isdigit():
		print('Error1: Number of non coding zones must be a number.')
		sys.exit(2)
	ref_file_idx = 21 + 2*nr_of_noncoding_reg	#index of reference file
	ref_file = info[ref_file_idx]		#name and address of reference file
	if not path.exists(ref_file):
		print('Error2: Cannot find reference file.')
		sys.exit(2)
	first_coordinate = info[ref_file_idx+1]	#first coordinate of rDNA unit
	if not str(first_coordinate).isdigit():
		print('Error3: First coordinate of rDNA unit must be a number.')
		sys.exit(2)
	last_coordinate = info[ref_file_idx+2]	#last_coordinate of rDNA unit
	if not str(last_coordinate).isdigit():
		print('Error4: Last coordinate of rDNA unit must be a number.')
		sys.exit(2)
	simulated_vcf_file = info[ref_file_idx+3]	#simulated vcf file
	simulated_rDNA_file = info[ref_file_idx+4]	#simulated rDNA file
	length_of_line = info[ref_file_idx+5]		#length of line of simulated rDNA fasta file
	if not str(length_of_line).isdigit():
		print('Error5: Length of line must be a number.')
		sys.exit(2)
	read_depth = info[ref_file_idx+6]		#coverage
	snp_error_rate = info[ref_file_idx+7]         #error rate of snps
	read_file_prefix = info[ref_file_idx+8]         #error rate of snps
	#bringing reference file ti linux standard
	process_cat = subprocess.Popen(['cat', ref_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	with open ('fixed.fasta', 'w') as fixedfile:
		subprocess.Popen(['tr', '-d', '\'\\r\''], stdin=process_cat.stdout, stdout=fixedfile, stderr=subprocess.PIPE).wait()
        #trim reference file
	subprocess.Popen([python, 'fasta_trimmer.py', 'fixed.fasta', 'ref_trm.fasta', '{0}'.format(first_coordinate), '{0}'.format(last_coordinate)]).wait()
	#generating the simulated vcf
	subprocess.Popen([python, 'variant_simulator.py', '{0}'.format(info)]).wait()
	#generating the simulated rDNA
	subprocess.Popen([python, 'genome_generator.py', 'ref_trm.fasta', '{0}'.format(simulated_vcf_file), '{0}'.format(simulated_rDNA_file), '{0}'.format(nr_of_repeats), '{0}'.format(length_of_line)]).wait()
	#activating the pIRS
	command = "export LD_LIBRARY_PATH=/usr/local/lib64:${LD_LIBRARY_PATH}"
	subprocess.Popen([command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	#generating reads
	devnull = open(os.devnull, 'w')	#to suppress screen output
	try:
		subprocess.Popen(['pirs', 'simulate', '{0}'.format(simulated_rDNA_file), '-x', '{0}'.format(read_depth), '-e', '{0}'.format(snp_error_rate), '-o', '{0}'.format(read_file_prefix)], stdout=devnull, stderr=devnull)
	except OSError:
		print('Error6: Cannot find pirs.')
		sys.exit(2)
	
#main function
def main():
	info = parser()
	variant_generator(info)

if __name__ == "__main__":
	python = 'python'+str(sys.version_info.major)
	config = sys.argv[1]                 #input reference file
	if not path.exists(config):
		print('Error0: Cannot find input file')
		sys.exit(2)
	main()


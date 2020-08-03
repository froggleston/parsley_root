import sys
import time
import subprocess
import os

#function to retrieve frequencies
def freq_retriever(chrom_pos):	#function to retrieve frequencies
	numparvar = len(chrom_pos)	#number of partial variants
	subprocess.Popen(['samtools', 'sort', 'bamfile.bam', '-o', 'bamfile_sorted.bam']).wait()
	subprocess.Popen(['samtools', 'index', 'bamfile_sorted.bam']).wait()
	for i in range(0, numparvar):	#creating array of commands to create vcf file
		process_call_variants = subprocess.Popen(['freebayes', '-f', 'ref_trm.fasta', '-F', '0.001', '--pooled-continuous', '--haplotype-length', '-1', '-r', chrom_pos[i], 'bamfile_sorted.bam'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		subprocess.Popen([python, 'header_remover.py', str(i), 'vcf_hr.vcf'], stdin=process_call_variants.stdout).wait()	#append on the file
        #Break multi-allelic sites
	process_breakmulti = subprocess.Popen(['vcfbreakmulti', 'vcf_hr.vcf'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #Decompose biallelic sites
	process_allelicprimitives = subprocess.Popen(['vcfallelicprimitives', '-kg'], stdin=process_breakmulti.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #Fix the field of alternate allelic depth
	process_fix_depth = subprocess.Popen(['./vcf_ad_fix.sh'], stdin=process_allelicprimitives.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #normalise the vcf
	process_normalise = subprocess.Popen(['bcftools', 'norm', '-f', 'ref_trm.fasta', '-m-'], stdin=process_fix_depth.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #relocate variants
	subprocess.Popen([python, 'variant_relocator.py', 'relocated_fr.vcf', str(shoulder_size+1), str(rDNA_unit_size+shoulder_size), '1'], stdin=process_normalise.stdout).wait()	
        #sort variants
	with open ('sorted_fr.vcf', 'w') as sortedfile:                      #file object to write on the file
		subprocess.Popen(['bcftools', 'sort', 'relocated_fr.vcf'], stdout=sortedfile, stderr=subprocess.PIPE).wait()
	subprocess.Popen([python, 'variant_unifier.py', 'sorted_fr.vcf', 'unified_fr.vcf', '1']).wait()
	subprocess.Popen([python, 'variant_unifier.py', 'unified_fr.vcf', 'unified_fr2.vcf', '0']).wait()      #unifying the relocated variants 
	with open ('filtered_fr.vcf', 'w') as filtered_fr_file:		#file object to write on the file	
		subprocess.Popen(['vcffilter', '-f', 'AO > 1 & AO / DP > 0.005', 'unified_fr2.vcf'], stdout=filtered_fr_file, stderr=subprocess.PIPE).wait()	#filtering the variants
	subprocess.Popen([python, 'data_updater.py', 'filtered.vcf', 'filtered_fr.vcf', vcf_outfile])   #updating the data

##function to parse vcf file
def parser_vcf_file():
	info = []
	with open('filtered.vcf', 'r') as file_handle:	#open vcf file for reading
		for line in file_handle:
			if line.startswith('#'):
				continue
			char_list = ''                          #initialising character list for reading
			AO_field = False                        #Its not time to fill AO field
			DP_field = False                        #its not time to fill DP field
			count = 0				#counts spaces in the line
			for index in range(len(line)):	#loop to parse each lin
				if count < 2:
					if not line[index].isspace():		#if parser doesn't encounter space
						char_list += line[index]	#add character to the character list
					else:					#if previusly value is inserted
						info.append(char_list)		#inserting value in info array
						char_list = ''			#initialising character list for reading
						count += 1
				elif line[index] == '=' and line[index-1] == 'O' and line[index-2] == 'A' and line[index-3] == ';':	#looking for pattern ;AO= 
					AO_field = True		#its time to fill AO field
				elif AO_field:
					if line[index].isdigit():
						char_list += line[index]        #add character to the character list
					else:
						info.append(char_list)          #inserting value in info array
						char_list = ''                  #initialising character list for reading
						AO_field = False		#AO field filled
				elif line[index] == '=' and line[index-1] == 'P' and line[index-2] == 'D' and line[index-3] == ';':     #looking for pattern ;DP=
					DP_field = True		#its time to fill DP field
				elif DP_field:
					if line[index].isdigit():
						char_list += line[index]        #add character to the character list
					else:
						info.append(char_list)          #inserting value in info array
						break
	return info

#main function
def main():
	info = parser_vcf_file()
	info_size = len(info)			#size of information array
	AO = []					#array to store number of reads representing alternate allele
	for i in range(2, info_size, 4):  	#loop to search AO field in info
		AO.append(int(info[i]))
	DP = []					#array to store total depth at variant loci
	for i in range(3, info_size, 4):  	#loop to search DP field
		DP.append(int(info[i]))
	numvariants = len(AO)			#number of vatiants taken from vcf file
	chrom_pos = []				#chromosome positions to search
	for i in range(numvariants):
		chrom_pos.append(int(info[4*i+1]) + shoulder_size - 1)	#lower bound
		chrom_pos.append(int(info[4*i+1]) + shoulder_size)	#upper bound
		if int(info[4*i+1]) < shoulder_size:                            #for left shoulder
			chrom_pos.append(int(info[4*i+1]) + rDNA_unit_size - 1)	#lower bound
			chrom_pos.append(int(info[4*i+1]) + rDNA_unit_size)	#upper bound
		elif rDNA_unit_size - int(info[4*i+1]) < shoulder_size:		#for right shoulder
			chrom_pos.append(max(1, shoulder_size + int(info[4*i+1]) - rDNA_unit_size - 1))	#lower bound
			chrom_pos.append(shoulder_size + int(info[4*i+1]) - rDNA_unit_size)		#upper bound
	chrom_pos.sort()		#sorting the array
	loop_break = True		#to start while loop
	while(loop_break):		#for filtering the close chromosome positions
		loop_break = False
		for i in range(2, len(chrom_pos), 2):
			if chrom_pos[i] <= chrom_pos[i-1] + 1:
				chrom_pos.pop(i-1)	#filter the upper bound
				chrom_pos.pop(i-1)	#filter the lower bound
				loop_break = True
				break
	chrom_rng = []			#for storing chrom positions as ranges
	for i in range(0, len(chrom_pos), 2):
		char_list = ''          #initiating charater list
		char_list += info[0]
		char_list += ':'
		char_list += str(chrom_pos[i])	#this is lower bound of locus
		char_list += '-'
		char_list += str(chrom_pos[i+1])  #this is upper bound of locus
		chrom_rng.append(char_list)				#storing in the array of chromosome positions
	freq_retriever(chrom_rng)	#retrieving the variant frequencies
			
if __name__ == "__main__":
	python = 'python'+str(sys.version_info.major)   #python version
	shoulder_size = int(sys.argv[1])	#shoulder size in fasta and bam file
	rDNA_unit_size = int(sys.argv[2])	#size of rDNA unit
	vcf_outfile = sys.argv[3]		#output vcf file
	main()

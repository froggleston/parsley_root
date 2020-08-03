"""Relocate Variants"""

import sys

#function to relocate variants in vcf file
def relocate_variants():
	VHI = 1000                                                                              #very high index in a vcf file line
	with open(rlcVCF, 'w') as file_handle2:                                      #open vcf file for writing
		for line in sys.stdin:                                                       #line loop for vcf file in reading
			if line.startswith('##contig'):
				(lineseg1, locus_nr_init_trm, locus_nr_fnl_trm, genome_length) = compute_coordinates_vcf(line) #computing coordinates of the trimmed vcf
				#writing line containing coordinates to trimmed vcf
				file_handle2.write(lineseg1 + ':' + str(locus_nr_init_trm) + '..' + str(locus_nr_fnl_trm) + ',length=' + str(genome_length) + '>' + '\n')
			elif not line.startswith('#'):                                                #if header is over
				chrName = ''
				for index in range(len(line)):
					if line[index] == ':':
						break
					chrName += line[index]
				(modstr, lineseg) = compute_coordinates_variants(VHI, line)					#computing coordinates of variants
				file_handle2.write(chrName + ':' + str(locus_nr_init_trm) + '..' + str(locus_nr_fnl_trm) + '\t' + modstr + lineseg) #print line with modiied coordinates
			else:
				file_handle2.write(line)                                                #exact copy of the line from the old file

#function to compute coordinates of the trimmed vcf
def compute_coordinates_vcf(line):
	lineseg1 = ''
	for index1 in range(0, len(line) - 1):      #loop to see from where initial locus number starts
		if line[index1] == ':':
			break
		lineseg1 += line[index1]
	locus1 = ''                                 #variable for storing initial locus number
	for index2 in range(index1 + 1, len(line) - 1):     #index from where initial locus number starts
		if line[index2].isdigit() == False:
			break
		locus1 += line[index2]                  #storing initial locus number
	locus_nr_init = int(locus1)                 #finding integer value
	genome_length = locus_end - locus_beg + 1               #computing total length of genome
	if mod_chr_coord == 1:
		locus_nr_init_trm = locus_nr_init + locus_beg - 1   #finding integer value of initial locus number of trimmed fasta file
		locus_nr_fnl_trm = locus_nr_init + locus_end - 1    #finding integer value of last locus number of trimmed fasta file
	else:
		locus_nr_init_trm = locus_nr_init			#keeping chromosome coordinate same for the trimmed fasta file
		locus_nr_fnl_trm = locus_nr_init + 2*(locus_beg-1) + genome_length - 1	#keeping chromosome coordinate same for the trimmed fasta file
	return (lineseg1, locus_nr_init_trm, locus_nr_fnl_trm, genome_length)

#function to compute coordinates of variants
def compute_coordinates_variants(VHI, line):
	for locus_index_in in range(0, VHI):                                    #loop to locate start index of DNA locus on the line in a vcf file
		if line[locus_index_in].isspace() == True:                          #DNA locus comes after space character
			break
	for locus_index_out in range(locus_index_in + 1, VHI):                  #loop to locate start index of DNA locus on the line in a vcf file
		if line[locus_index_out].isspace() == True:                         #to check line locus where digits end
			break
	orgstr = ''                                                             #initialize string of an original number
	for locus in range(locus_index_in + 1, locus_index_out):                #loop to construct number string
		orgstr += line[locus]
	locus_nr = int(orgstr)                                                  #convert string into a number
	if (locus_nr < locus_beg):
		modstr = str(locus_end - 2*locus_beg + locus_nr + 2)
	elif (locus_nr <= locus_end):
		modstr = str(locus_nr - locus_beg + 1)                              #reset the coodinate i.e. locus number in the modified string
	else: modstr = str(locus_nr - locus_end)
	lineseg = ''                                                       #third line segment of rest of the line after the string of locus number
	for locus in range(locus_index_out, len(line)):
		lineseg += line[locus]
	return(modstr, lineseg)

#main function
def main():
	relocate_variants()

if __name__ == "__main__":
	locus_beg = int(sys.argv[2])
	locus_end = int(sys.argv[3])
	rlcVCF = sys.argv[1]
	mod_chr_coord = int(sys.argv[4])	#modify choromosome coordinates (1/0) (True/False)
	main()


"Program to generate genome of strain"

import sys
from random import randint

#function to parse the fasta file
def parser_REF(REF):
	DNA_seq = []
	with open(REF, 'r') as file_handle:                     #open reference fasta file for reading
		for line in file_handle:                       #loop of the lines in the ref fasta file
			if line.startswith('>'):                        #condition to identify first line of fasta
				lineseg1 = ''
				for index1 in range(1, len(line) - 1):      #loop to store chromosome name by skipping starting character
					if (line[index1].isspace()):
						if ':' in lineseg1:
                                                	break
						else:
							lineseg1 = ''
					else:
						lineseg1 += line[index1]
				lineseg2 = ''                                               #storing last segment
				for index2 in range (index1, len(line) - 1):
					lineseg2 += line[index2]
			else:
				if line[len(line)-1] == "\n":                   #carraige return is additional locus of the line
					for index in range(len(line)-1):
						DNA_seq.append(line[index])
				else:                                   #condition to identify last line of fasta
					for index in range(len(line)):
						DNA_seq.append(line[index])
		GS = len(DNA_seq)		#length of one unit rDNA(Genome size)
		for i in range(CN-1):
			for index in range(GS):
				DNA_seq.append(DNA_seq[index])	#total length of rDNA
	return (GS, DNA_seq, lineseg1, lineseg2)

#function to parse the VCF file
def parser_VCF(VCF, lineseg1):
	coordinate = []                                 #variant coordinates
	CNV = []                                        #copy number variants
	REF_list = []                                   #array of reference bases
	ALT_list = []                                   #array of alt bases
	with open(VCF, 'r') as file_handle:             #open VCF file for reading
		for line in file_handle:                        #loop of the lines in the VCF file
			if line.startswith(lineseg1):           #condition to identify line showing variant of VCF
				digit_list = ''
				REF_seq = ''
				digit = False			#digit is not parsed yet
				alpha = False			#alphabetic character is not parsed yet
				point = False			#point is not parsed yet
				for index1 in range(len(lineseg1), len(line) - 1):
					if line[index1].isdigit():
						digit_list += line[index1]
						digit = True			#digit parsed
					elif line[index1].isalpha():		#to check if its alphabetical character
						REF_seq += line[index1]
						alpha = True			#alphabetic character inserted
					elif digit and alpha:		#if alphabetic and numeric characters and point are parsed
                                		break
				coordinate.append(int(digit_list))
				REF_list.append(REF_seq)
				digit_list = ''
				ALT_seq = ''
				digit = False				#any digit is not captured yet
				alpha = False				#any aphabetic character is not captured yet
				for index2 in range(index1, len(line) - 1):
					if line[index2].isdigit():
						digit_list += line[index2]
						digit = True		#digit captured
					elif line[index2].isalpha():
						ALT_seq += line[index2]
						alpha = True		#alphabetic character captured
					elif digit and alpha:		#if both alphabetic and numeric characters are captured
                                		break
				CNV.append(int(digit_list))
				ALT_list.append(ALT_seq)
	return (coordinate, CNV, REF_list, ALT_list)

#function to generate whole DNA of chromosome of the variant strain
def strain_DNA_generator(DNA_seq, coord_glob, ref_list_exp, alt_list_exp):
	NV = len(coord_glob)	#total number of variants in all repeats
	numsnps = 0
	numdels = 0
	numins = 0
	dl = False
	ins = False
	for i in range(NV-1, -1, -1):
		if len(ref_list_exp[i]) == len(alt_list_exp[i]):
			DNA_seq[coord_glob[i]-1] = alt_list_exp[i]	#this is SNP/pSNP
			numsnps += 1
		elif len(ref_list_exp[i]) > len(alt_list_exp[i]):
			del DNA_seq[coord_glob[i]:coord_glob[i]+len(ref_list_exp[i])-1]	#this is deletion
			numdels += 1
		else:
			for j in range(len(alt_list_exp[i])-1, 0, -1):
				DNA_seq.insert(coord_glob[i], alt_list_exp[i][j])	#this is insertion
			numins += 1
			if len(alt_list_exp[i]) > 10:
				seq = []
				k = -1
				for j in range(coord_glob[i]-1, coord_glob[i]+33):
					seq.append('')
					k += 1
					seq[k] = DNA_seq[j]
	return (DNA_seq)

#function for assigning rDNA_units to variants
def assigning_rDNA_units(CN, GS, coordinate, CNV, REF_list, ALT_list):
	NV = len(coordinate)            		##number of variants
	base_muted = [False for i in range(GS*(CN+1))]	##base muted	
	CNV_array = []
	for j in range(NV):
		CNV_list = []
		m = 0
		while m < CNV[j]:
			if CNV[j] == CN:
				i = m
			else:
				i = randint(0, CN-1)
			cont = False
			if len(REF_list[j]) == len(ALT_list[j]):        #if this is a snp
                                if base_muted[coordinate[j]+i*GS]:      #check if base is already muted
                                        cont = True
			elif len(REF_list[j]) > len(ALT_list[j]):       #if this is deletion
                                for k in range(1, len(REF_list[j])):
                                        if base_muted[coordinate[j]+i*GS+k]:    #check if bases are already muted
                                                cont = True
                                                break
			else:                                           #this is insertion
                                if base_muted[coordinate[j]+i*GS+1]:    #check if bases are already muted
                                        cont = True
			if cont:                        #skip if bases are muted
				continue 
			CNV_list.append(str(i))
			if len(REF_list[j]) == len(ALT_list[j]):                 ##this is snp
				base_muted[coordinate[j]+i*GS] = True
			elif len(REF_list[j]) > len(ALT_list[j]):                       #if this is deletion
				for k in range(1, len(REF_list[j])):                 #mute the bases
					base_muted[coordinate[j]+i*GS+k] = True
			else:                                                   #if this is insertion and snp
				base_muted[coordinate[j]+i*GS+1] = True
			m += 1
		CNV_array.append(CNV_list)
	return CNV_array

#function for expanding vcf
def vcf_expander(CN, GS, coordinate, CNV_array, REF_list, ALT_list):
	NV = len(coordinate)	#number of variants
	ref_list_exp = []	#new expanded reference list
	alt_list_exp = []	#new expanded alt list
	coord_glob = []		#global coordinate spreading over all rDNA repeats
	total_variants = 0	#total number of variants
	for i in range(CN):
		for j in range(NV):
			if str(i) in CNV_array[j]:
				ref_list_exp.append('')
				alt_list_exp.append('')
				coord_glob.append(coordinate[j]+i*GS)
				ref_list_exp[total_variants] = REF_list[j]
				alt_list_exp[total_variants] = ALT_list[j]
				total_variants += 1
	return (coord_glob, ref_list_exp, alt_list_exp)

#function for writing genome file
def file_writing(GEN, DNA_seq, lineseg1, lineseg2, NBL):
	with open(GEN, 'w') as file_handle:      #open vcf file for writing
		file_handle.write('>' + lineseg1 + lineseg2) #header of genome file
		size = len(DNA_seq)					#size of genome
		for index in range(size):
			if index%NBL == 0:
				file_handle.write('\n')
			file_handle.write(DNA_seq[index])     
                
#main function
def main():
	(GS, DNA_seq, lineseg1, lineseg2) = parser_REF(REF)	#function to parse reference sequence
	(coordinate, CNV, REF_list, ALT_list) = parser_VCF(VCF, lineseg1)		#function to parse VCF file
	CNV_array = assigning_rDNA_units(CN, GS, coordinate, CNV, REF_list, ALT_list)	#function to assign rDNA units to variants
	(coord_glob, ref_list_exp, alt_list_exp) = vcf_expander(CN, GS, coordinate, CNV_array, REF_list, ALT_list)	#function to expand vcf to all rDNA repeats
	DNA_seq = strain_DNA_generator(DNA_seq, coord_glob, ref_list_exp, alt_list_exp)	#function to generate whole DNA (all repeats) of chromosome of the variant strain
	file_writing(GEN, DNA_seq, lineseg1, lineseg2, NBL)


if __name__ == "__main__":
	REF = sys.argv[1]               #input reference file
	VCF = sys.argv[2]               #iutput vcf file
	GEN = sys.argv[3]        	#output genome file
	CN = int(sys.argv[4])           #copy numbers
	NBL = int(sys.argv[5])		#number of bases in one line
	main()


"""Unify Variants"""

import sys
from  math import log10

#function to unify variants in vcf file
def unify_variants():
	with open(orgVCF, 'r') as file_handle1:                                          #open vcf file for reading
		with open(uniVCF, 'w') as file_handle2:                                      #open vcf file for writing
			first_line = True				#first line to start
			for line in file_handle1:
				if line.startswith('#'):		#if not header
					file_handle2.write(line)
				else:
					#extracting variant information
					if first_line:
						(lineseg11, pos1, lineseg21, ref1, alt1, qlt1, lineseg31, ao1, lineseg41, dp1, lineseg51, ro1, lineseg61, lineseg71) = extract_info(line)
						first_line = False		#first line over
					else:
						(lineseg12, pos2, lineseg22, ref2, alt2, qlt2, lineseg32, ao2, lineseg42, dp2, lineseg52, ro2, lineseg62, lineseg72) = extract_info(line)
						if (pos1 == pos2 and ref1 == ref2 and alt1 == alt2 and dec_var == 1 and dp1 == dp2): #if two variants originated from compound variant match
							ao1 += ao2                      				#add number of reads of alternative alleles of two variants
							q = (10**(-qlt1/10)*dp1 + 10**(-qlt2/10)*dp2)/(dp1 + dp2)
							if q == 0:
								qlt1 = min(qlt1, qlt2)
							else:
								qlt1 = -10*log10(q) 	#compute quality of merged variant
						elif (pos1 == pos2 and ref1 == ref2 and alt1 == alt2 and dec_var == 0 and dp1 != dp2):		#if two variants match because of relocation
							ao1 += ao2				#add number of reads of alternative alleles of two variants
							ro1 += ro2                              #add number of reads of reference alleles of two variants
							dp1 += dp2                              #add coverage of two variants
							q = (10**(-qlt1/10)*dp1 + 10**(-qlt2/10)*dp2)/(dp1 + dp2)
							if q == 0:
								qlt1 = min(qlt1, qlt2)
							else:
								qlt1 = -10*log10(q)     #compute quality of merged variant
						else:								#if two variants do not match
							file_handle2.write(lineseg11 + pos1 + lineseg21 + ref1 + '\t' + alt1 + '\t' + str(qlt1) + lineseg31 + str(ao1) + lineseg41 + str(dp1) + ';DPB=' + str(dp1) + lineseg51 + str(ro1) + lineseg61 + str(dp1) + ':' + str(dp1) + ',' + str(ao1) + ':' + str(ro1) + ':' + lineseg71 + '\n')	#print the last line
							#transfering information of line2 to line1
							lineseg11 = lineseg12
							pos1 = pos2
							lineseg21 = lineseg22
							ref1 = ref2
							alt1 = alt2
							qlt1 = qlt2
							lineseg31 = lineseg32
							ao1 = ao2
							lineseg41 = lineseg42
							dp1 = dp2
							lineseg51 = lineseg52
							ro1 = ro2
							lineseg61 = lineseg62
							lineseg71 = lineseg72
			file_handle2.write(lineseg11 + pos1 + lineseg21 + ref1 + '\t' + alt1 + '\t' + str(qlt1) + lineseg31 + str(ao1) + lineseg41 + str(dp1) + ';DPB=' + str(dp1) + lineseg51 + str(ro1) + lineseg61 + str(dp1) + ':' + str(dp1) + ',' + str(ao1) + ':' + str(ro1) + ':' + lineseg71 + '\n')   #print the last line

#function to extract variant information
def extract_info(line):
	lineseg1 = ''
	for index1 in range(0, len(line) - 1):      #loop to skip chromosome name
		lineseg1 += line[index1]
		if line[index1].isspace() == True:
			break
	pos = ''						#variable to save position (locus number) of variant
	lineseg2 = ''
	for index2 in range(index1 + 1, len(line) - 1):     #index from where locus number of variant  starts
		if line[index2].isdigit() == False:
			lineseg2 += line[index1]
			break
		pos += line[index2]                  #storing locus number of variant
	for index1 in range(index2+1, len(line) - 1):      #loop to skip id field
		lineseg2 += line[index1]
		if line[index1].isspace() == True:
			break
	ref = ''                                           #variable to save reference allele of variant
	for index2 in range(index1 + 1, len(line) - 1):     #index from where reference allele of variant  starts
		if line[index2].isspace() == True:
			break;
		ref += line[index2]                  #storing reference allele of variant
	alt = ''					#variable to save alternate allele of variant
	for index1 in range(index2 + 1, len(line) - 1):      #index from where alternate allele of variant starts
		if line[index1].isspace() == True:
			break
		alt += line[index1]
	qlt = ''						#quality value
	lineseg3 = ''						#variable to save quality of variant
	for index2 in range(index1 + 1, len(line) - 1):		#index from where quality of variant  starts
		if line[index2].isspace() == True:
			lineseg3 += line[index1]
			break
		qlt += line[index2]                  #storing quality of variant
	for index1 in range(index2 + 1, len(line) - 1):      #loop to skip filter field
		lineseg3 += line[index1]
		if line[index1].isspace() == True:
			break
	for index2 in range(index1 + 1, len(line) - 1):     #loop to skip info field up to AO field
		lineseg3 += line[index2]
		if line[index2]  == '=' and line[index2-1]  == 'O' and line[index2-2]  == 'A':
			break
	ao = ''
	lineseg4 = ''
	for index1 in range(index2 + 1, len(line) - 1):     #loop to to store AO  field
		if line[index1].isdigit() == False:
			lineseg4 += line[index1]
			break
		ao += line[index1]
	for index2 in range(index1 + 1, len(line) - 1):     #loop to skip info field up to DP field
		lineseg4 += line[index2]
		if line[index2]  == '=' and line[index2-1]  == 'P' and line[index2-2]  == 'D':
			break
	dp = ''
	for index1 in range(index2 + 1, len(line) - 1):     #loop to store DP  field
		if line[index1].isdigit() == False:
			break
		dp += line[index1]
	lineseg5 = ''
	for index2 in range(index1 + 1, len(line) - 1):     #loop to skip info field up to DPB field
		if line[index2] == ';':
			lineseg5 += line[index2]
			break
	for index1 in range(index2 + 1, len(line) - 1):      #loop to skip info field up to RO field
		lineseg5 += line[index1]
		if line[index1]  == '=' and line[index1-1]  == 'O' and line[index1-2]  == 'R' and line[index1-3]  == ';':
			break
	ro = ''
	lineseg6 = ''
	for index2 in range(index1 + 1, len(line) - 1):     #loop to store RO  field
		if line[index2].isdigit() == False:
			lineseg6 += line[index2]
			break
		ro += line[index2]
	for index1 in range(index2 + 1, len(line) - 1):      #loop to skip rest of the info field
		lineseg6 += line[index1]
		if line[index1].isspace()  == True:
			break
	for index2 in range(index1 + 1, len(line) - 1):      #loop to skip format field
		lineseg6 += line[index2]
		if line[index2].isspace() == True:
			break
	for index1 in range(index2 + 1, len(line) - 1):      #loop to skip GT field
		lineseg6 += line[index1]
		if line[index1] == ':':
			break
	for index2 in range(index1 + 1, len(line) - 1):      #loop to skip DP field
		if line[index2] == ':':
			break
	for index1 in range(index2 + 1, len(line) - 1):      #loop to skip AD field
		if line[index1] == ':':
			break
	for index2 in range(index1 + 1, len(line) - 1):      #loop to skip AO field
		if line[index2] == ':':
			break
	for index1 in range(index2 + 1, len(line) - 1):      #loop to skip RO field
		if line[index1] == ':':
			break
	lineseg7 = ''
	for index2 in range(index1 + 1, len(line) - 1):      #loop to skip dp field
		lineseg7 += line[index2]
	return (lineseg1, pos, lineseg2, ref, alt, float(qlt), lineseg3, int(ao), lineseg4, int(dp), lineseg5, int(ro), lineseg6, lineseg7)

def main():
	unify_variants()

if __name__ == "__main__":
	orgVCF = sys.argv[1]		#original/input vcf
	uniVCF = sys.argv[2]		##unified/output vcf
	dec_var = int(sys.argv[3])	#decision variable to decide whether to unify relocated variants or decomposed variants
	main()

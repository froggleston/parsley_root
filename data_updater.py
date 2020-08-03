"""data updater"""

import sys

#function to unify variants in vcf file
def update_data():
	with open(outVCF, 'w') as writing_file:				#open file for writing
		with open(orgVCF, 'r') as original_file:			#open file for reading
			for line in original_file:
				if line.startswith('#'):
					writing_file.write(line)
					continue
				(lineseg_org1, pos_org, lineseg_org2, ref_org, alt_org, qlt_org, lineseg_org3, ao_org, lineseg_org4, dp_org, lineseg_org5, ro_org, lineseg_org6, lineseg_org7, lineseg_org8) = extract_info(line)
				(success, lineseg_cor1, pos_cor, lineseg_cor2, ref_cor, alt_cor, qlt_cor, lineseg_cor3, ao_cor, lineseg_cor4, dp_cor, lineseg_cor5, ro_cor, lineseg_cor6, lineseg_cor7, lineseg_cor8) = match(pos_org, ref_org, alt_org)
				if success:
					writing_file.write(lineseg_cor1 + pos_cor + lineseg_cor2 + ref_cor + '\t' + alt_cor + '\t' + str(qlt_cor) + '\t' + lineseg_cor3 + str(ao_cor) + lineseg_cor4 + str(dp_cor) + ';DPB=' + str(dp_cor) + lineseg_cor5 + str(ro_cor) + lineseg_cor6 + str(dp_cor) + ':' + str(dp_cor) + ',' + str(ao_cor) + ':' + str(ro_cor) + ':' + lineseg_cor7 + str(ao_cor) + lineseg_cor8 + '\n')
				else:
					writing_file.write(line)

#to match variants
def match(pos_org, ref_org, alt_org):
	success = False
	with open(corVCF, 'r') as corrected_file:        #open file for reading
		for line in corrected_file:
			if not line.startswith('#'):            #if header
				(lineseg_cor1, pos_cor, lineseg_cor2, ref_cor, alt_cor, qlt_cor, lineseg_cor3, ao_cor, lineseg_cor4, dp_cor, lineseg_cor5, ro_cor, lineseg_cor6, lineseg_cor7, lineseg_cor8) = extract_info(line)
				if(pos_cor == pos_org and ref_cor == ref_org and alt_cor == alt_org):
					success = True
					break
	return (success, lineseg_cor1, pos_cor, lineseg_cor2, ref_cor, alt_cor, qlt_cor, lineseg_cor3, ao_cor, lineseg_cor4, dp_cor, lineseg_cor5, ro_cor, lineseg_cor6, lineseg_cor7, lineseg_cor8)

#function to extract variant information
def extract_info(line):
	lineseg1 = ''
	for index1 in range(0, len(line) - 1):      #loop to skip chromosome name
		lineseg1 += line[index1]
		if line[index1].isspace():
			break
	pos = ''                                   #variable to save position (locus number) of variant
	lineseg2 = ''
	for index2 in range(index1 + 1, len(line) - 1):     #index from where locus number of variant  starts
		if not line[index2].isdigit():
			lineseg2 += line[index1]
			break
		pos += line[index2]                  #storing locus number of variant
	for index1 in range(index2+1, len(line) - 1):      #loop to skip id field
		lineseg2 += line[index1]
		if line[index1].isspace():
			break
	ref = ''                                           #variable to save reference allele of variant
	for index2 in range(index1 + 1, len(line) - 1):     #index from where reference allele of variant  starts
		if line[index2].isspace():
			break;
		ref += line[index2]                  #storing reference allele of variant
	alt = ''                                    #variable to save alternate allele of variant
	for index1 in range(index2 + 1, len(line) - 1):      #index from where alternate allele of variant starts
		if line[index1].isspace():
			break
		alt += line[index1]
	qlt = ''                                       #variable to save quality of variant
	for index2 in range(index1 + 1, len(line) - 1):     #index from where quality of variant  starts
		if line[index2].isspace():
			break
		qlt += line[index2]                  #storing quality of variant
	lineseg3 = ''				
	for index1 in range(index2 + 1, len(line) - 1):      #loop to store filter field
		lineseg3 += line[index1]
		if line[index1].isspace():
			break
	for index2 in range(index1 + 1, len(line) - 1):     #loop to skip info field up to AO field
        	lineseg3 += line[index2]
        	if line[index2]  == '=' and line[index2-1]  == 'O' and line[index2-2]  == 'A':
            		break
	ao = ''
	for index1 in range(index2 + 1, len(line) - 1):     #loop to store AO  field
		if not line[index1].isdigit():
			break
		ao += line[index1]
	lineseg4 = ';'
	for index2 in range(index1 + 1, len(line) - 1):     #loop to store info field up to DP field
		lineseg4 += line[index2]
		if line[index2]  == '=' and line[index2-1]  == 'P' and line[index2-2]  == 'D':
			break
	dp = ''
	for index1 in range(index2 + 1, len(line) - 1):     #loop to store DP  field
		if not line[index1].isdigit():
			break
		dp += line[index1]
	for index2 in range(index1 + 1, len(line) - 1):     #loop to skip DPB  field
		if line[index2] == ';':
			break
	lineseg5 = ';'
	for index1 in range(index2 + 1, len(line) - 1):      #loop to store info field up to RO field
		lineseg5 += line[index1]
		if line[index1]  == '=' and line[index1-1]  == 'O' and line[index1-2]  == 'R' and line[index1-3]  == ';':
			break
	ro = ''
	for index2 in range(index1 + 1, len(line) - 1):     #loop to store RO  field
		if not line[index2].isdigit():
			break
		ro += line[index2]
	lineseg6 = ';'
	count = 0
	for index1 in range(index2 + 1, len(line) - 1):      #loop to store rest of the info field and upto GT format field
		lineseg6 += line[index1]
		if line[index1]  == ':':
			count += 1
			if count == 8:
				break
	for index2 in range(index1 + 1, len(line) - 1):      #loop to skip DP, AD, RO formate fields
		if line[index2] == ':':
			count += 1
		if count == 11:
			break;
	lineseg7 = ''
	for index1 in range(index2 + 1, len(line) - 1):      #loop to store QR field
		lineseg7 += line[index1]
		if line[index1] == ':':
			break
	for index2 in range(index1 + 1, len(line) - 1):      #loop to skip AO field
		if line[index2] == ':':
			break
	lineseg8 = ':'	
	for index1 in range(index2 + 1, len(line) - 1):      #loop to store QA & GL field
		lineseg8 += line[index1]
	return (lineseg1, pos, lineseg2, ref, alt, float(qlt), lineseg3, int(ao), lineseg4, int(dp), lineseg5, int(ro), lineseg6, lineseg7, lineseg8)

def main():
        update_data()

if __name__ == "__main__":
	orgVCF = sys.argv[1]		#original/input vcf
	corVCF = sys.argv[2]		#vcf with corrected frequencies
	outVCF =  sys.argv[3]		#output vcf
	main()


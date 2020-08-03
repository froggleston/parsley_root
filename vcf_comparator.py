#comparing vcf files

import sys

def parser(VCF, f):
	with open(VCF, 'r') as file_handle:        #open file for readin
		info = []			#complete information array
                locus = []                      #coordinates of all loci
                ref = []                        #reference allele array
                alt = []                        #alternate allele array
		if f == 's':			#if simulation file
                	cnv = []                        #copy numbers of this variant (array)
		else:				#if prediction file
               		AO = []                        #Number of alternate allele
                	DP = []                         #depth of reads
	        for line in file_handle:
			space_done =  [False for i in range(8)]   #space not arrived yet
                	locus_temp = ''                 #coordinate of one locus
                	ref_temp = ''                   #one reference allele
                	alt_temp = ''                   #one alternate allele
			if f == 's':			#if simulation file
                		cnv_temp = ''                   #copy numbers of one variant
			else:				#if prediction file
                        	AO_temp = ''                   #frequency of one variant
                        	DP_temp = ''
        	        if line.startswith('#'):            #if header
				continue
			for index in range(len(line)):
	                	if not space_done[0]:
					if line[index].isspace():
						space_done[0] = True		#first space arrived
				elif not space_done[1]:
					if line[index].isdigit():
						locus_temp += line[index]
					else:
						space_done[1] = True		#second space arrived
				elif not space_done[2]:
                                        if line[index].isspace():
                                                space_done[2] = True            #third space arrived
				elif not space_done[3]:
					if line[index].isalpha():
						ref_temp += line[index]
					else:
						space_done[3] = True		#fourth space arrived
				elif not space_done[4]:
                                        if line[index].isalpha():
                                                alt_temp += line[index]
                                        else:
                                                space_done[4] = True            #fifth space arrived
				elif f == 's':					#if simulation file
					if line[index].isdigit():
                                		cnv_temp += line[index]

					else:
						locus.append(locus_temp)
						ref.append(ref_temp)
						alt.append(alt_temp)
						cnv.append(cnv_temp)
						break
                                elif not space_done[5]:
                                        if line[index] == '=' and line[index-1] == 'O' and line[index-2] == 'A' and line[index-3] == ';':
                                                space_done[5] = True
                                elif not space_done[6]:
                                        if line[index].isdigit():
                                                AO_temp += line[index]
                                        else:
                                                space_done[6] = True
                                elif not space_done[7]:
                                        if line[index] == '=' and line[index-1] == 'P' and line[index-2] == 'D' and line[index-3] == ';':
                                                space_done[7] = True
                                elif line[index].isdigit():
                                        DP_temp += line[index]
                                else:
                                        locus.append(locus_temp)
                                        ref.append(ref_temp)
                                        alt.append(alt_temp)
                                        AO.append(AO_temp)
                                        DP.append(DP_temp)
                                        break
		info.append(locus)
		info.append(ref)
		info.append(alt)
		if f == 's':				#if simulation file
			info.append(cnv)
		else:					#if prediction file
			info.append(AO)
                        info.append(DP)
	return info

#comparing two vcfs
def VCF_comparer(sim_info, prd_info, cn):
	numvariants_sim = len(sim_info[0])
	numvariants_prd = len(prd_info[0])
	allele_done =  [False for i in range(numvariants_prd)]  
	comp_info = []					#comparative vcf info
	for i in range(numvariants_sim):		#variants from simulation vcf
		for j in range(3):
			comp_info.append(sim_info[j][i])
		comp_info.append(str(format(float(sim_info[3][i])/float(cn), '.3f')))	#frequency in simulated vcf
		match = False
		for j in range(numvariants_prd):		#matching variants from predictive vcf
			if prd_info[0][j] == sim_info[0][i] and  prd_info[1][j] == sim_info[1][i] and  prd_info[2][j] == sim_info[2][i]:
				comp_info.append(str(format(float(prd_info[3][j])/float(prd_info[4][j]), '.3f'))) #frequency in predictive vcf
				comp_info.append('Exact Match')
				allele_done[j] = True
				match = True
				break
			elif int(prd_info[0][j]) >= int(sim_info[0][i]) - 5 and int(prd_info[0][j]) <= int(sim_info[0][i]) and len(sim_info[1][i]) != len(sim_info[2][i]) and len(prd_info[1][j]) - len(prd_info[2][j]) == len(sim_info[1][i]) - len(sim_info[2][i]):
				comp_info.append(str(format(float(prd_info[3][j])/float(prd_info[4][j]), '.3f')))
				comp_info.append('Left Aligned')
				allele_done[j] = True
				match = True
				break
		if not match:
			comp_info.append('0.000')	#frequency in predictive vcf of non matching variants
			comp_info.append('Not Found')
        for i in range(numvariants_prd):		#variants from predictive vcf which did not match
		if allele_done[i]:			#if already matched
			continue
                for j in range(3):
                        comp_info.append(prd_info[j][i])
		comp_info.append('0.000')			#this did not match with the simulated vcf
                comp_info.append(str(format(float(prd_info[3][j])/float(prd_info[4][j]), '.3f'))) #frequency in predictive vcf
		comp_info.append('False Positive')
		
	return comp_info	

def main():
	sim_info = parser(simVCF, 's')
	prd_info = parser(prdVCF, 'p')
	comp_info = VCF_comparer(sim_info, prd_info, cn)
	with open(compVCF, 'w') as file_handle:
		file_handle.write('Locus'.ljust(25) + 'Ref'.ljust(25) + 'alt'.ljust(25)  + 'frq_sim'.ljust(25) + 'frq_prd'.ljust(25) + 'Comment' + '\n')
		for i in range(len(comp_info)):
			if (i+1)/6 < float(i+1)/float(6):
				file_handle.write(comp_info[i].ljust(25))
			else:
				file_handle.write(comp_info[i] + '\n')


if __name__ == "__main__":
        simVCF = sys.argv[1]	#simulated vcf
	cn = float(sys.argv[2])	#copy numbers
        prdVCF = sys.argv[3]	#predicted vcf
        compVCF =  sys.argv[4]	#comparative vcf
        main()


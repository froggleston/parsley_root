"Program to generate variations"

import sys
import math
import random
from random import randint
from random import choice
import numpy as np
				 
#function to parse the fasta file
def parser():
	DNA_seq = ''
	with open('ref_trm.fasta', 'r') as file_handle:			#open reference fasta file for reading
		for line in file_handle:			#loop of the lines in the ref file
			if line.startswith('>'):		#condition to identify first line of fasta
				lineseg1 = ''
				for index1 in range(1, len(line) - 1):		#loop to store chromosome name by skipping starting character 
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
				if line[len(line)-1] == "\n":			#carraige return is additional locus of the line
					for index in range(len(line)-1):
						DNA_seq += line[index]
				else:					#condition to identify last line of fasta
					for index in range(len(line)):
                                        	DNA_seq += line[index]
	return (DNA_seq, lineseg1, lineseg2)

#function to generate snps
def SNP_generator(DNA_seq, genome_size, var_status, times_base_muted):
	vz = random.uniform(0, 1)	#variant zone i.e. to know whether it lies in non-coding zone
	if vz <= prob_of_var_in_noncoding_zone:
		(locus, var_size, fail) = coord_selector_noncoding_zone(var_status, times_base_muted, 0, [1, 1], 'snp')
	else:
		(locus, var_size, fail) = coord_selector_rDNA(var_status, times_base_muted, 0, [1, 1], genome_size, 'snp')
	mutant = mutant_selector(DNA_seq, locus)
	frq = frequency_assigner(var_status, nr_of_copies-times_base_muted[locus])
	return (locus, mutant, frq, fail)

#coordinate selector for deletion in non-homopolymeric region and non-coding zone
def coord_selector_noncoding_zone(var_status, times_base_muted, times_base_muted_indel, var_range, var_type):
	loop = True
	fail = False				#coordinate selection is not failed yet
	itr = 0					#number of iterations
	while(loop):
		itr += 1			#increment in number of iterations
		loop = False
		ncz_id = randint(0, nr_of_noncoding_reg-1)  #selection of non-coding zone
		var_size = randint(var_range[0], var_range[1])    #random selection of variation size for deletion
		locus = randint(LB_for_noncoding_reg[ncz_id], UB_for_noncoding_reg[ncz_id]-var_size)   #locus to be muted
		if var_type != 'snp' and (times_base_muted_indel[locus-1] > 0 or times_base_muted_indel[locus+var_size] > 0):   #to avoid continuous indel
			loop = True
		for j in range(var_size):
			if (times_base_muted[locus+j] >  part_var_freq_range[1] or (var_status == 'fixed' and times_base_muted[locus+j] > 0)):
				loop = True
				break
		if (itr > 1000):		#if more than 1000 iterations
			fail = True		#coordinate selector failed
			break
	return (locus, var_size, fail)

#coordinate selector for deletion and insertion for any zone (coding or non-coding)
def coord_selector_rDNA(var_status, times_base_muted, times_base_muted_indel, var_range, genome_size, var_type):
	loop = True
	fail = False                            #coordinate selection is not failed yet
	itr = 0                                 #number of iterations
	while(loop):
		itr += 1                        #increment in number of iterations
		loop = False
		var_size = randint(var_range[0], var_range[1])    #random selection of variation size
		locus = randint(1, genome_size-var_size)   #locus to be muted
		if var_type != 'snp' and (times_base_muted_indel[locus-1] > 0 or times_base_muted_indel[locus+var_size] > 0):	#to avoid continuous indel
			loop = True
		for j in range(var_size):
			if (times_base_muted[locus+j] > part_var_freq_range[1] or (var_status == 'fixed' and times_base_muted[locus+j] > 0)):
				loop = True
				break
		if (itr > 1000):                #if more than 1000 iterations
			fail = True             #coordinate selector failed
			break
	return (locus, var_size, fail)

#this function selects the mutant base
def mutant_selector(DNA_seq, locus):
	bases = ['A', 'C', 'G', 'T']    #DNA bases
	while(True):                    #loop for selection of mutant
		mutant = choice(bases)  #selecting a mutant
		if mutant != DNA_seq[locus]:    #if mutant is different from the reference
			break
	return mutant

#function to generate deletions
def DEL_generator(DNA_seq, genome_size, var_status, times_base_muted, times_base_muted_del, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones):
	variant_zone = random.uniform(0, 1)       #variant zone i.e. to know whether it is in non-coding zone
	if variant_zone <= prob_of_var_in_noncoding_zone:
		variant_zone = random.uniform(0, 1)       #to see whether variant is in homopolymeric zone
		if variant_zone <= prob_homo_indel:
			#coordinate selector in homopolymeric zone
			(locus, var_size, fail) = coord_selector_homo_zone(var_status, times_base_muted, times_base_muted_del, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, 3, 'del')
		else:
			#coordinate selector in noncoding zone
			(locus, var_size, fail) = coord_selector_noncoding_zone(var_status, times_base_muted, times_base_muted_del, [min_del_size, max_del_size], 'del')
	else:
		#coordinate selector in full rDNA
		(locus, var_size, fail) = coord_selector_rDNA(var_status, times_base_muted, times_base_muted_del, [min_del_size, max_del_size], genome_size, 'del')
	freq = frequency_assigner(var_status, nr_of_copies-times_base_muted[locus])
	return (locus, var_size, freq, fail)

#function to generate insertions
def INS_generator(DNA_seq, genome_size, var_status, times_base_muted, times_base_muted_ins, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, Largest_homo_size):
	bases = ['A', 'C', 'G', 'T']
	inserts = ''
	vz = random.uniform(0, 1)       #variant zone i.e. to know whether it is in non-coding zone
	if vz <= prob_of_var_in_noncoding_zone:
		vz = random.uniform(0, 1)	#to see variant is in homopolymeric region
		if vz <= prob_homo_indel:
			#coordinate selector in homopolymeric zone
			(locus, var_size, fail) = coord_selector_homo_zone(var_status, times_base_muted, times_base_muted_ins, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, Largest_homo_size, 'ins')
			for j in range(var_size):
				inserts += DNA_seq[locus]
		else:
			#coordinate selector in noncoding zone
			(locus, var_size, fail) = coord_selector_noncoding_zone(var_status, times_base_muted, times_base_muted_ins, [min_ins_size, max_ins_size], 'ins')
			for j in range(var_size):
				inserts += choice(bases)
	else:
		#coordinate selector in full rDNA unit
		(locus, var_size, fail) = coord_selector_rDNA(var_status, times_base_muted, times_base_muted_ins,[min_ins_size, max_ins_size], genome_size, 'ins')
		for j in range(var_size):
			inserts += choice(bases)
	frq = frequency_assigner(var_status, nr_of_copies-times_base_muted[locus])
	return (locus, inserts, frq, fail)

#generate coordinate of homopolymeric zone
def coord_selector_homo_zone(var_status, times_base_muted, times_base_muted_indel, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, bound, var_type):
	##to compute slots available in homopolymeric regions/zones
	nslots = compute_slots(LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, times_base_muted_indel, times_base_muted)
	var_size = 0				#initializing size of indel
	itr = 0					#number of iterations
	fail = False				#coordinate selector is not failed yet
	loop = True
	while(loop):
		itr += 1			#increment in number of iterations
		if var_status == 'partial':			##if partial indel
			##setting probability distribution for likelihood selection of homopolymeric zone
			prb_dist = []			#probability distribution
			if nslots[0] == 0:		##if no slots avaiable
				prb_dist.append(0)
			else:
				prb_dist.append(nr_of_copies-nslots[0]+1)	##probability range for the first homopolymeric zone
			for i in range(1,  nr_of_homo_zones):				##probability distribution for the rest of homopolymeric zones
				if nslots[i] == 0:	#if less vacant slots than minimum allowable frequency in this homopolymeric zone
					prb_dist.append(prb_dist[i-1]+0)
				else:
					prb_dist.append(prb_dist[i-1]+nr_of_copies-nslots[i]+1)
			rnd_prb = randint(1, prb_dist[nr_of_homo_zones-1])	##random selection from probability distribution
			for i in range(nr_of_homo_zones):			##selection of homopolymeric zone based on probability distribution
				if rnd_prb <= prb_dist[i]:	
					homo_zone = i
					break
			homo_size = UB_for_homo_zone[homo_zone] - LB_for_homo_zone[homo_zone] + 1                 #computing size of homopolymer
			itr2 = 0
			while(True):
				itr2 += 1
				var_size = indel_size(var_type, homo_size, bound) 	##computing size of indel
				if (times_base_muted[LB_for_homo_zone[homo_zone] + var_size] == times_base_muted[LB_for_homo_zone[homo_zone] + var_size - 1]):	#to stop same sized variation
					break
				if itr2 > 10:
					break
			if itr2 <= 10:
				break
		else:
			while True:
				itr += 1
				homo_zone = randint(0,  nr_of_homo_zones-1)
				if nslots[homo_zone] == nr_of_copies: ##if all slots are vacant
					homo_size = UB_for_homo_zone[homo_zone] - LB_for_homo_zone[homo_zone] + 1                 #computing size of homopolymer
					var_size = indel_size(var_type, homo_size, bound)	##computing size of indel
					loop = False
					break
				if itr > 100:
					break
		if (itr > 100):		#if more than 100 iterations
			fail = True		#coordinate selector failed
			break
	return (LB_for_homo_zone[homo_zone], var_size, fail)

##to compute slots available in homopolymeric regions/zones
def compute_slots(LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, times_base_muted_indel, times_base_muted):
	nslots = []             ##number of available slots
	for i in range(nr_of_homo_zones):
		if times_base_muted_indel[LB_for_homo_zone[i]-1] > 0:
			nslots.append(0)        ##no slot available if there is indel prior to start of the zone
			continue
		nfslots = 0      #number of filled slots
		for j in range(LB_for_homo_zone[i], UB_for_homo_zone[i]+1):
			nfslots = max(nfslots, times_base_muted[j])
		if nr_of_copies-nfslots < part_var_freq_range[0]:	#if number of available slots are less than the allowable lower bound on frequency of the partial variation
			nslots.append(0)
		else:
			nslots.append(nr_of_copies-nfslots)
	return nslots

##computing size of indel
def indel_size(var_type, homo_size, bound):
	if (var_type == 'del'):
		var_size = randint(1, homo_size - bound)  #computing size of deletion
	else:
		var_size = randint(1, bound - homo_size)  #computing size of insertion
	return var_size

#function to create vcf file
def file_writing(DNA_seq, lineseg1, lineseg2, SNP_list, SNP_coord, SNP_freq, DEL_coord, DEL_size, DEL_freq, INS_array, INS_coord, INS_freq):
	SNP_done =  [False for j in range(len(SNP_coord))]   #SNP not chosen yet
	DEL_done =  [False for j in range(len(DEL_coord))]   #DEL not chosen yet
	INS_done =  [False for j in range(len(INS_coord))]   #INS not chosen yet
	nVars = len(SNP_coord) + len(DEL_coord) + len(INS_coord)
	with open(VCF, 'w') as file_handle:      #open VCF file for writing
		file_handle.write('##Simulated Variants for ' + lineseg1 + lineseg2 + '\n') #header of VCF file
		file_handle.write('#CHROM' + '\t' + 'POS' + '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'FREQ' + '\t' + 'TYPE' + '\n')
		SNP = 0					#no SNP has been written
		DEL = 0					#no DEL has been written
		INS = 0					#no INS has been written
		for j in range(nVars):
			coordinate = 10000
			for i in range(len(SNP_coord)):
				if SNP_coord[i] < coordinate and not SNP_done[i]:
					coordinate = SNP_coord[i]
					coord = i
					Type = 'SNP'
			for i in range(len(DEL_coord)):
				if DEL_coord[i]-1 < coordinate and not DEL_done[i]:
					coordinate = DEL_coord[i]-1
					coord = i
					Type = 'DEL'
			for i in range(len(INS_coord)):
				if INS_coord[i]-1 < coordinate and not INS_done[i]:
					coordinate = INS_coord[i]-1
					coord = i
					Type = 'INS'
			if Type == 'SNP':
				SNP_done[coord] = True
			elif Type == 'DEL':
				DEL_done[coord] = True
			else:
				INS_done[coord] = True
			file_handle.write(lineseg1 + '\t' + str(coordinate+1) + '\t' + '.' + '\t' + DNA_seq[coordinate]) #fasta coordinate starts from 1 python coordinate starts from 0
			if Type == 'SNP':
				file_handle.write('\t' + SNP_list[coord])
				FREQ = SNP_freq[coord]
			elif Type == 'DEL':
				for i in range(1, DEL_size[coord]+1):
					file_handle.write(DNA_seq[coordinate + i])
				file_handle.write('\t' + DNA_seq[coordinate])
				FREQ = DEL_freq[coord]

			else:								#this is for insertion
				file_handle.write('\t' + DNA_seq[coordinate])	
				file_handle.write(INS_array[coord])
				FREQ = INS_freq[coord]
			file_handle.write('\t' + str(FREQ))
			if FREQ == nr_of_copies:
				file_handle.write('\t' + Type + '\n')
			else:
				file_handle.write('\t' + 'p' + Type + '\n')

#function todecide frequency of the variation
def frequency_assigner(var_status, frq_available):
	if var_status == 'fixed' or frq_available <= part_var_freq_range[0]:           #This should cover all copy numbers
		frq = frq_available	
	else:                           #This should cover part of copy numbers
		#random frequency based on beta distributionn and limits set by user
		frq = max(min(int(math.ceil(np.random.beta(alpha, beta, size=None)*frq_available)), part_var_freq_range[1]), part_var_freq_range[0])
	return frq

#function to locate homopolymeric zones
def find_homopolymeric_zones(DNA_seq):
	LB_for_homo_zone = []	#lower bound homopolymeric zone
	UB_for_homo_zone = [] 	#upper bound homopolymeric zone
	largest_homo_size = 0		#size of largest homopolymer
	for j in range(nr_of_noncoding_reg):
		T = 0           #number of consecutive Ts
		A = 0           #number of consecutive As
		homo = False    #homopolymeric zone is not started yet
		for i in range (LB_for_noncoding_reg[j]-1, UB_for_noncoding_reg[j]):
			T = compute_homopolymeric_base(DNA_seq[i], 'T', T)
			A = compute_homopolymeric_base(DNA_seq[i], 'A', A)
			if T == min_homo_size or A == min_homo_size:
				LB_for_homo_zone.append(i-min_homo_size+1)	#minimum size for homopolymeric zone attained
				homo = True	#we are in a homopolymeric zone
			if homo and T < min_homo_size and A < min_homo_size:
				UB_for_homo_zone.append(i-1)	#homopolymeric zone ended 1 base before
				homo = False	#end of homopolymeric tract
				largest_homo_size = max(largest_homo_size, UB_for_homo_zone[len(UB_for_homo_zone)-1] - LB_for_homo_zone[len(LB_for_homo_zone)-1] + 1)
	return (LB_for_homo_zone, UB_for_homo_zone, largest_homo_size)

#function to compute homopolymeric bases
def compute_homopolymeric_base(DNA, base, n):
	if DNA == base:
		n += 1
	else:
		n = 0
	return n

#generation of variants
def variant_generator(DNA_seq, lineseg1, lineseg2, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, Largest_homo_size):
	genome_size = len(DNA_seq)
	times_base_muted = [0 for i in range(genome_size-min_ins_size+max_ins_size)] #size increased for variation which could go beyond standard size due to insertion
	times_base_muted_del = [0 for i in range(genome_size-min_ins_size+max_ins_size)]     #exclusively for deletions
	times_base_muted_ins = [0 for i in range(genome_size-min_ins_size+max_ins_size)]     #exclusively for insertions
	SNP_coord = []                                   #coordinates of SNP
	DEL_coord = []                                   #coordinates of deletion
	INS_coord = []                                   #coordinates of insertion
	SNP_list = []                                   #list of SNP mutant bases
	DEL_size = []                                   #list of deletion sizes
	INS_list = []                                   #list of inserted bases
	SNP_frq = []                                    #list of frequency of SNP
	DEL_frq = []                                    #list of frequency of deletion
	INS_frq = []                                    #list of frequency of insertion
	fail = False
	i = 0
	while i < nr_of_var:                           ##loop to generate all variants
                s_vr = random.uniform(0, 1)             #variation selector
                if s_vr <= cum_prob_SNP:
                        (locus,  mutant, frq, fail) = SNP_generator(DNA_seq, genome_size, 'fixed', times_base_muted) #function for the generation of SNPs
                elif s_vr <= cum_prob_part_SNP:
                        (locus, mutant, frq, fail) = SNP_generator(DNA_seq, genome_size, 'partial', times_base_muted) #function for the generation of pSNPs
                elif s_vr <= cum_prob_DEL:      #function for generation DELs
                        (locus, size, frq, fail) = DEL_generator(DNA_seq, genome_size, 'fixed', times_base_muted, times_base_muted_del, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones)
                elif s_vr <= cum_prob_part_DEL:     #function for generation pDELs
                        (locus, size, frq, fail) = DEL_generator(DNA_seq, genome_size, 'partial', times_base_muted, times_base_muted_del, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones)
                elif s_vr <= cum_prob_INS:      #function for the generation of INS
                        (locus, inserts, frq, fail) = INS_generator(DNA_seq, genome_size, 'fixed', times_base_muted, times_base_muted_ins, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, Largest_homo_size)
                else:                   #function for the generation of pINS
                        (locus, inserts, frq, fail) = INS_generator(DNA_seq, genome_size, 'partial', times_base_muted, times_base_muted_ins, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, Largest_homo_size)
                if (fail):
                        continue
                if s_vr <= cum_prob_part_SNP:
                        SNP_coord.append(locus)          #appending to SNP coordinate list
                        SNP_list.append(mutant)         #appending to mutant list
                        SNP_frq.append(frq)             #appending to SNP frequency
                        times_base_muted[locus] += frq        #adding number of copies of variations to this locus
                elif s_vr <= cum_prob_part_DEL:
                        DEL_coord.append(locus)          #appending to DEL coordinate list
                        DEL_size.append(size)           #appending to deletion size
                        DEL_frq.append(frq)             #appending to deletion frequency
                        for j in range(size):
                                times_base_muted[locus+j] += frq              #adding number of copies of variations to these loci
                                times_base_muted_del[locus+j] += frq          #for deletion only
                else:
                        INS_coord.append(locus)          #appending to INS coordinate list
                        INS_list.append(inserts)        #appending to insertion list
                        INS_frq.append(frq)             #appending to insertion frequency
                        size = len(inserts)
                        for j in range(size):
                                times_base_muted[locus+j] += frq      #adding number of copies of variations to this locus
                                times_base_muted_ins[locus+j] += frq  #for insertion only
                i += 1
	file_writing(DNA_seq, lineseg1, lineseg2, SNP_list, SNP_coord, SNP_frq, DEL_coord, DEL_size, DEL_frq, INS_list, INS_coord, INS_frq)       #file writing

#main function
def main():
	(DNA_seq, lineseg1, lineseg2) = parser()    #parsing the reference genome file
	(LB_for_homo_zone, UB_for_homo_zone, largest_homo_size) = find_homopolymeric_zones(DNA_seq)
	nr_of_homo_zones = len(LB_for_homo_zone)					#number of homopolymeric regions
	Largest_homo_size = max(max_homo_size, largest_homo_size + 1)			#largest allowable homopolymer size in this program
	variant_generator(DNA_seq, lineseg1, lineseg2, LB_for_homo_zone, UB_for_homo_zone, nr_of_homo_zones, Largest_homo_size)

class take_command_line_args():
	def __init__(self, info):
		global nr_of_copies, nr_of_noncoding_reg, LB_for_noncoding_reg, UB_for_noncoding_reg, NB, nr_of_var, min_del_size, max_del_size, min_ins_size, max_ins_size, min_homo_size, max_homo_size, part_var_freq_range, cum_prob_SNP, cum_prob_part_SNP, cum_prob_DEL, cum_prob_part_DEL, cum_prob_INS, cum_prob_INS, prob_of_var_in_noncoding_zone, prob_homo_indel, alpha, beta, VCF
		nr_of_copies = int(info[0])                     #copy numbers or number of repeats of rDNA array
		nr_of_noncoding_reg = int(info[1])                               #number of non-coding regions
		LB_for_noncoding_reg = []                                                #lower bounds of non-coding regios
		for i in range(2, 2 + nr_of_noncoding_reg):                     #lower bound array
			LB_for_noncoding_reg.append(int(info[i]))                 #add lower bound	
		UB_for_noncoding_reg = []                                         #upper bounds of non-coding regios
		NB = 2*nr_of_noncoding_reg                                       #total number of bounds (upper & lower)
		for i in range(2 + nr_of_noncoding_reg, 2 + NB):                        #upper bound array
			UB_for_noncoding_reg.append(int(info[i]))                 #add upper bound
		nr_of_var = int(info[2 + NB])                          #total number of variants
		min_del_size = int(info[3 + NB])                    #lower bound on deletion size
		max_del_size = int(info[4 + NB])                    #upper bound on deletion size
		min_ins_size = int(info[5 + NB])                    #lower bound on insertion size
		max_ins_size = int(info[6 + NB])                    #upper bound on insertion size
		min_homo_size = int(info[7 + NB])                       #lower bound on homopolymer size
		max_homo_size = int(info[8 + NB])                    #upper bound on homopolymer size
		part_var_freq_range = []                                       #partial variation frequency range
		part_var_freq_range.append(int(info[9 + NB]))                  #lower bound on partial variation frequency range
		part_var_freq_range.append(int(info[10 + NB]))                 #upper bound on partial variation frequency range
		cum_prob_SNP = float(info[11 + NB])                     #commulative probability of SNP
		cum_prob_part_SNP = cum_prob_SNP + float(info[12 + NB])             #commulative probability of pSNP
		cum_prob_DEL = cum_prob_part_SNP + float(info[13 + NB])             #commulative probability of deletion
		cum_prob_part_DEL =  cum_prob_DEL + float(info[14 + NB])             #commulative probability of partial deletion
		cum_prob_INS = cum_prob_part_DEL + float(info[15 + NB])             #commulative probability of insertion
		cum_prob_part_INS = cum_prob_INS +  float(info[16 + NB])            #commulative probability of  partial insertion
		prob_of_var_in_noncoding_zone = float(info[17 + NB])                    #probability of variation in non-coding zone
		prob_homo_indel = float(info[18 + NB])                #probability of homopolymeric indels
		alpha = float(info[19 + NB])                    #parameter of beta distribution
		beta = float(info[20 + NB])                     #parameter of beta distribution
		VCF =  info[24 + NB]                            #file writing
	def execute(self):
		main()
	

if __name__ == "__main__":
	info = sys.argv[1].split(', ')		#split string into list of elements
	info[0] = info[0].replace('[','',1)	#remove bracket from first character
	info[len(info)-6] = info[len(info)-6].replace('\'','',1)	#remove open quote from character
	info[len(info)-6] = info[len(info)-6].replace('\'','',1)	#remove close quote from character
	take_command_line_args(info).execute()

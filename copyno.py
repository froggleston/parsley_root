#program to estimate number of rDNA copies

import sys
from scipy import stats
import numpy as np
from sklearn.linear_model import LinearRegression

#set up sliding window function
def slideFunct(data, window, step):
	total = len(data)
	spots = []
	for i in range (0, total-window, step):
		spots.append(i)
	result = []
	for i in range (len(spots)):
		result.append(np.mean(data[spots[i]:spots[i]+window]))
	return(result)

#parameter file parser
def parser_prm(REF):
	info = []
	with open(REF, 'r') as file_handle:                     #open parameter file for reading
		for line in file_handle:                        #loop of the lines in the config file
			read = False
			char_list = ''
			if (not line.startswith('#') and not len(line) < 2):
				if (line[len(line)-1] == "\n"):
					line_len = len(line) - 1
				else:
					line_len = len(line)
				for index in range(line_len):          #loop to parse each lin
					if (line[index] == '='):
						read = True
					elif (read):
						if (not line[index].isspace()):
							char_list += line[index]        #add value to the character list
				info.append(char_list)
	return info

#depth file parser
def parser_dfile(REF):
	info = []
	with open(REF, 'r') as file_handle:                     #open reference fasta file for reading
		for line in file_handle:                        #loop of the lines in the config file
			col = 0
			char_list = ''
			for index in range(len(line)):                  #loop to parse each line
				if col < 2:				#skip two colums
					if (line[index].isspace()):     #if space or tab
						col += 1
					continue
				elif (line[index].isdigit()):		#if digit
					char_list += line[index]        #add value to the character list
				else:					#if not digit
					info.append(char_list)          #inserting value in info
	return info

##examine the read depth across the whole of Chromosome XII
def examine_whole_chromosome(REF, rdnastart, rdnaend, bufer, window, step):
	wgsdata12 = parser_dfile(REF)		#parsing the depth file
	wgsdata12 = np.array(wgsdata12).astype(np.int)
	startdata12 = []
	for j in range(rdnastart-bufer):
		startdata12.append(wgsdata12[j])
	enddata12 = []
	for j in range(rdnaend+bufer, len(wgsdata12)):
		enddata12.append(wgsdata12[j])
	notrdna12 = np.concatenate((startdata12, enddata12))
	sw12 = slideFunct(notrdna12, window, step)
	length = len(sw12)
	sw12index = []
	for j in range(length):
		sw12index.append(j)
	slope, intercept, r_value, p_value, std_err = stats.linregress(sw12index,sw12)
	return intercept

#examine the read depth across the rDNA single unit, forcing the regression line to be horizontal
def examine_rdna_unit(REF, bufer, window, step):
	rdnadata = parser_dfile(REF)
	rdnadata = np.array(rdnadata).astype(np.int)
	for j in range(bufer):
		rdnadata[len(rdnadata)-2*bufer+j] += rdnadata[j]              #adding right shoulder
		rdnadata[bufer+j] += rdnadata[len(rdnadata)-bufer+j]        #adding left shoulder
		rdnadata2 = []
		for j in range(bufer+1, len(rdnadata)-bufer):
			rdnadata2.append(rdnadata[j])
		sw12 = slideFunct(rdnadata2, window, step)
	return sw12

def examine(whole_genome_depth_profile, rdnastart, rdnaend, buffer1, rDNA_depth_profile, shoulder, window, step):
	#examine the read depth across the whole of Chromosome XII
	intercept = examine_whole_chromosome(whole_genome_depth_profile, rdnastart, rdnaend, buffer1, window, step)
	#examine the read depth across the rDNA single unit, forcing the regression line to be horizontal
	sw12 = examine_rdna_unit(rDNA_depth_profile, shoulder, window, step)
	return int(round(np.mean(sw12)/intercept))

#main function
def main():
  return examine(whole_genome_depth_profile, rdnastart, rdnaend, buffer1, rDNA_depth_profile, shoulder, window, step)

if __name__ == "__main__":
	whole_genome_depth_profile = sys.argv[1]	#whole genome depth profile text file
	rdnastart = int(sys.argv[2])				#rDNA start locus in whole genome
	rdnaend =  int(sys.argv[3])				#rDNA end locus in whole genome
	buffer1 = int(sys.argv[4])				#buffer between rDNA and non rDNA region
	rDNA_depth_profile = sys.argv[5]		#depth profile of rDNA region only
	shoulder = int(sys.argv[6])				#its equal to shoulder size
	window = int(sys.argv[7])				#window size
	step = int(sys.argv[8])				#step size
	main()

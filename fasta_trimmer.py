"""Trim a Fasta file."""
import sys

#function to trim fasta file
def trim_fasta(orgFasta, trmFasta, locus_beg, locus_end):
	info = []
	(lenLine, numLines_beg, numLines_end, rem_beg, rem_end) = parameters(orgFasta, locus_beg, locus_end)	#computing parameter
	with open(orgFasta, 'r') as file_handle1:          		    #open original fasta file for reading
		with open(trmFasta, 'w') as file_handle2:      		    #open trimmed fasta file for writing
			count_lines = 0                                     #counting lines in original file
			count_seq = 0                                       #count sequence number of locus within line
			lineseg = ''                                        #initializing the line segment of trimmed fasta file
			for line in file_handle1:                           #loop of the lines in the original file
				if line.startswith('>'):                    #condition to identify first line of fasta
					(lineseg1, lineseg2, locus_nr_init_trm, locus_nr_fnl_trm) = compute_coordinates(line)	#computing coordinates of the trimmed fasta
					file_handle2.write(lineseg1 + ':' + str(locus_nr_init_trm) + '..' + str(locus_nr_fnl_trm) + lineseg2 + '\n')	#writing line containing coordinates to trimmed fasta file
					info.append(lineseg1 + ':' + str(locus_nr_init_trm) + '..' + str(locus_nr_fnl_trm) + lineseg2 + '\n')	#for file stream
				else:
					count_lines += 1
					if count_lines > numLines_end:         	#if number of lines exceeds end limit of trimmed file
						break
					if count_lines < numLines_beg:            	#if number of lines is less than the start limit of the trimmed file
						continue
					if(count_lines == numLines_beg)&(count_lines == numLines_end):
						for locus in range(rem_beg, rem_end + 1): 	#loop to write starting line of DNA sequence to trimmed file
							(lineseg, count_seq, info) = file_writing(file_handle2, line, lineseg, locus, count_seq, lenLine, info)
					elif count_lines == numLines_beg:        #if number of lines is equal to where start locus of trimmed file present
						for locus in range(rem_beg, lenLine): #loop to write starting line of DNA sequence to trimmed file
							(lineseg, count_seq, info) = file_writing(file_handle2, line, lineseg, locus, count_seq, lenLine, info)
					elif count_lines == numLines_end:      #if number of lines is equal to where end locus of trimmed file present
						for locus in range(0, rem_end + 1):         # loop to write end line of DNA sequence to trimmed file
							(lineseg, count_seq, info) = file_writing(file_handle2, line, lineseg, locus, count_seq, lenLine, info)
					else:
						for locus in range(0, lenLine):         #loop to write middle lines of DNA sequence to trimmed file
							(lineseg, count_seq, info) = file_writing(file_handle2, line, lineseg, locus, count_seq, lenLine, info)
			if(count_seq < lenLine):                            #writing last line of DNA sequence if it is less than full length
				file_handle2.write(lineseg)
				info.append(lineseg)	#for file stream
	return info

#function to compute coordinates of the trimmed fasta file
def compute_coordinates(line):
	lineseg1 = ''
	for index1 in range(0, len(line) - 1):      #loop to see from where initial locus number starts
        	if line[index1] == ':':
            		break
        	lineseg1 += line[index1]
	locus1 = ''                                 #variable for storing initial locus number
	for index2 in range(index1 + 1, len(line) - 1):	#index from where initial locus number starts
		if line[index2].isdigit() == False:
			break
		locus1 += line[index2]                  #storing initial locus number
	locus_nr_init = int(locus1)                 #finding integer value
	locus_nr_init_trm = locus_nr_init + locus_beg - 1   #finding integer value of initial locus number of trimmed fasta file
	locus_nr_fnl_trm = locus_nr_init + locus_end - 1	#finding integer value of last locus number of trimmed fasta file
	for index1 in range(index2 + 1, len(line) - 1):
		if line[index1].isspace() == True:
			break
	lineseg2 = ''						#storing last segment
	for index2 in range (index1, len(line) - 1):
		lineseg2 += line[index2]
	return (lineseg1, lineseg2, locus_nr_init_trm, locus_nr_fnl_trm)

#function to compute parameters of the DNA sequence
def parameters(orgFasta, locus_beg, locus_end):
	with open(orgFasta, 'r') as file_handle1:      #open original fasta file for reading
		lenLine = lengthOfLine(file_handle1)       #length of a line in an original fasta file
		numLines_beg = int(locus_beg/lenLine) + 1  #line number in original file where beginner locus of trimmed file exists
		numLines_end = int(locus_end/lenLine) + 1  #line number in original file where end locus of trimmed file exist
		rem_beg = locus_beg%lenLine - 1            #sequence number of beginner locus within line of original file
		if rem_beg < 0:
			rem_beg += lenLine                     #remainder 0 means last locus
			numLines_beg -= 1			   #line is not broken therefore it equals integer division
		rem_end = locus_end%lenLine - 1            #sequence number of end locus within line of original file
		if rem_end < 0:
			rem_end += lenLine                     #remainder 0 means last locus
			numLines_end -= 1			   #line is not broken therefore itequals integer division
		file_handle1.close()                       #close the original fasta file
	return(lenLine, numLines_beg, numLines_end, rem_beg, rem_end)

#function to write DNA sequence in a trimmed file
def file_writing(file_handle2, line, lineseg, locus, count_seq, lenLine, info):
	lineseg += line[locus]                              #adding DNA locus to line segment
	count_seq += 1                                      #incrementing sequence number counter of the line segment
	if(count_seq == lenLine):                           #if counter achieves the full length of line
		file_handle2.write(lineseg+'\n')                #write on the trimmed file
		info.append(lineseg+'\n')			#for file stream
		lineseg = ''                                    #initilize the next line segment of trimmed file
		count_seq = 0                                   #initialize the sequence counter of the next line in the trimmed file
	return (lineseg, count_seq, info)

#computing length of the line in the original file
def lengthOfLine(file_handle1):
	for line in file_handle1:                       #line loop in the original file
		if line.startswith('>'):                    #skip first line of original file
			continue
		return len(line) - 1

def main():
	trim_fasta(orgFasta, trmFasta, locus_beg, locus_end)
  
if __name__ == "__main__":
	locus_beg = int(sys.argv[3])	#begining locus
	locus_end = int(sys.argv[4])	#end locus
	orgFasta = sys.argv[1]		#original fasta/reference sequence file
	trmFasta = sys.argv[2]		#trimed fasta/output file
	main()

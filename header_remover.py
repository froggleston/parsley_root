#this is a program to remove extra headers from vcf file
import sys
#import io
#import subprocess

##function to parse vcf file
def parser_vcf_file(f):
	info = ''
	for line in sys.stdin:
		if f == 'w':				#if this is first process
			for index in range(len(line)):  #loop to parse each line
				info += line[index]
		elif not line.startswith('#'):		#record only nonheader lines
			for index in range(len(line)): 
				info += line[index]
	return info

#main function
def main():
	if(i == 0):
		f = 'w'
	else:
		f = 'a'
	info = parser_vcf_file(f)
	with open(vcf_outfile, f) as file_handle:      	#open vcf file for writing
		file_handle.write(info)

if __name__ == "__main__":
	vcf_outfile = sys.argv[2]               	#output vcf file
	i = int(sys.argv[1])				#sequence number of the function called
	main()


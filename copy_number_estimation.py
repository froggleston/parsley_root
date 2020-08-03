import sys, os, subprocess
import configparser

from itertools import islice

from fasta_trimmer import trim_fasta
from copyno import examine

def parser(parameters):
  config = configparser.ConfigParser()
  config.read(parameters)
  return config

def max_line_length(infile):
  max_length = 0
  max_line_length = ''
  with open(infile, 'r') as inf:
    while True:
      fastq_line_block = list(islice(inf, 4))
      if not fastq_line_block:
        break

      if(len(fastq_line_block[1]) > max_length):
        max_length = len(fastq_line_block[1])
        max_line_length = fastq_line_block[1]
  return max_length

def copy_number_estimator(config):
  cne_ref_params = config['cne_ref_params']
	whole_genome_ref = cne_ref_params['Name and address of whole genome fasta file']		#whole genome fasta file
  rDNA_start_locus = int(cne_ref_params['rDNA start locus number'])		#start locus of rDNA in whole genomw fasta file
  rDNA_end_locus = int(cne_ref_params['rDNA end locus number'])		#end locus of rDNA in whole genomw fasta file
  rDNA_buffer = int(cne_ref_params['buffer between rDNA and non-rDNA region'])		#buffer between rDNA and non rDNA

  cne_rdna_params = config['cne_rdna_params']
  rDNA_first_coord = int(cne_rdna_params['First Coordinate of rDNA unit'])		#first coordinate of rDNA in rDNA fasta file
  rDNA_last_coord = int(cne_rdna_params['Last Coordinate of rDNA unit'])		#last coordinate of rDNA in rDNA fasta file
  rDNA_ref = info[cne_rdna_params['Name and address of reference fasta file']]			#fasta file of rDNA only

  cne_read_params = config['cne_read_params']
  fwd_reads_file = cne_read_params['forward_read_ref']			#fastq file of forward reads
  rvs_reads_file = cne_read_params['reverse_read_ref']			#fastq file of reverse reads

  cne_slide_params = config['cne_slide_params']
  window_size = int(cne_slide_params['window size'])		#window size for computing average read depth
  step_size = int(cne_slide_params['step size'])		#step size for computing average read depth

	#change files to linux format
  with open('rDNA_fixed.fasta', 'w') as rDNA_fixedfile:
    with open(rDNA_ref) as rdna_ref_file:
      for line in rdna_ref_file:
        line.replace('\r\n', '\n')
        rDNA_fixedfile.write(line)

  with open('whole_genome_fixed.fasta', 'w') as whole_genome_fixedfile:
    with open(whole_genome_ref) as wg_ref_file:
      for line in wg_ref_file:
        line.replace('\r\n', '\n')
				whole_genome_fixedfile.write(line)

	#find maximum read length
  max_read_length_fwd = max_line_length(fwd_reads_file)   #maximum read length in forward read file
  max_read_length_rvs = max_line_length(rvs_reads_file)   #maximum read length in reverse read file

  max_read_length = max(max_read_length_fwd, max_read_length_rvs) #greater of forward and reverse reads
  shoulder_size = max_read_length - 1
  rDNA_first_coord = rDNA_first_coord - shoulder_size
  rDNA_last_coord = rDNA_last_coord + shoulder_size

  #trim the rDNA fasta file
  # subprocess.Popen(['python', 'fasta_trimmer.py', 'rDNA_fixed.fasta', 'rDNA_trm.fasta', '{0}'.format(rDNA_first_coord), '{0}'.format(rDNA_last_coord)]).wait()
  trim_fasta('rDNA_fixed.fasta', 'rDNA_trm.fasta', rDNA_first_coord, rDNA_last_coord)

  #index reference files
  subprocess.Popen(['bwa', 'index', 'rDNA_trm.fasta'], stderr=subprocess.PIPE)
  subprocess.Popen(['bwa', 'index', 'whole_genome_fixed.fasta'], stderr=subprocess.PIPE).wait()

  #map the reads
  devnull = open(os.devnull, 'w')
  with open ('rDNA_samfile.sam', 'w') as rDNA_samfile:
    subprocess.Popen(['bwa', 'mem', 'rDNA_trm.fasta', fwd_reads_file, rvs_reads_file], stdout=rDNA_samfile, stderr=devnull)	#to suppress screen output
	with open ('whole_genome_samfile.sam', 'w') as whole_genome_samfile:
    subprocess.Popen(['bwa', 'mem', 'whole_genome_fixed.fasta', fwd_reads_file, rvs_reads_file], stdout=whole_genome_samfile, stderr=devnull).wait()

  #sam to bam
  with open ('rDNA_bamfile.bam', 'w') as rDNA_bamfile:
    subprocess.Popen(['samtools', 'view', '-bS', 'rDNA_samfile.sam'], stdout=rDNA_bamfile, stderr=subprocess.PIPE)
  with open ('whole_genome_bamfile.bam', 'w') as whole_genome_bamfile:
    subprocess.Popen(['samtools', 'view', '-bS', 'whole_genome_samfile.sam'], stdout=whole_genome_bamfile, stderr=subprocess.PIPE).wait()

  #sort bam
  with open ('rDNA_bamsort.bam', 'w') as rDNA_bamsort:
    subprocess.Popen(['samtools', 'sort', 'rDNA_bamfile.bam'], stdout=rDNA_bamsort, stderr=subprocess.PIPE)
  with open ('whole_genome_bamsort.bam', 'w') as whole_genome_bamsort:
    subprocess.Popen(['samtools', 'sort', 'whole_genome_bamfile.bam'], stdout=whole_genome_bamsort, stderr=subprocess.PIPE).wait()

  #compute depth profile
  with open ('rDNA_depth.txt', 'w') as rDNA_depthfile:
    subprocess.Popen(['samtools', 'depth', '-aa', 'rDNA_bamsort.bam'], stdout=rDNA_depthfile, stderr=subprocess.PIPE)
  with open ('whole_genome_depth.txt', 'w') as whole_genome_depthfile:
    subprocess.Popen(['samtools', 'depth', '-aa', 'whole_genome_bamsort.bam'], stdout=whole_genome_depthfile, stderr=subprocess.PIPE).wait()

	#copy number estimate
	#copy_number_estimate = subprocess.Popen(['python3.6', 'copyno.py', 'whole_genome_depth.txt', '{0}'.format(rDNA_start_locus), '{0}'.format(rDNA_end_locus), '{0}'.format(rDNA_buffer), 'rDNA_depth.txt', '{0}'.format(shoulder_size), '{0}'.format(window_size), '{0}'.format(step_size)], stdout=subprocess.PIPE)
	#copy_number_estimate = copy_number_estimate.stdout.readlines()
  copy_number_estimate = examine('whole_genome_depth.txt', rDNA_start_locus, rDNA_end_locus, rDNA_buffer, 'rDNA_depth.txt', shoulder_size, window_size, step_size)
  return copy_number_estimate

#main function
def main(parameters):
  cne = copy_number_estimator(parser(parameters))
  print(cne)

if __name__ == "__main__":
  parameters = sys.argv[1]                #input parameter file
  main(parameters)

#!/usr/bin/python

import sys
import os
import getopt

def main():
	params = parseArgs()
	
	#build dictionary of reference contigs
	reference=dict()
	skip=0
	read=0
	print("Reading reference genome",params.ref)
	for contig in read_fasta(params.ref):
		if len(contig[1]) >= params.minlen:
			reference[contig[0]] = contig[1]
			read+=1
		else:
			skip+=1
			continue
	print("Done!\nTotal contigs read:",str(read))
	if skip > 0:
		print("Contigs skipped below min length:",str(skip))
	
	#read in coverage for sample 1
	

#Read samples as FASTA. Generator function
def read_fasta(fas):
	"""generator function loops through multifasa 
	records, returning tuples in the form of 
	[FastaHeader, FastaSequence]"""
	if os.path.exists(fas):
		with open(fas, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					#print(line)
					if line[0] == ">": #Found a header line
						#If we already loaded a contig, yield that contig and
						#start loading a new one
						if contig:
							yield([contig,seq]) #yield
							contig = "" #reset contig and seq
							seq = ""
						split_line = line.split()
						contig = (split_line[0].replace(">",""))
					else:
						seq += line
				#Iyield last sequence, if it has both a header and sequence
				if contig and seq:
					yield([contig,seq])
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)	

#Object to parse command-line arguments
class parseArgs():
	"""
	Class provides text command-line menu and holds parameters
	needed for the run
	"""
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'h1:2:r:o:znd:m:', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.ref=None
		self.sam1=None
		self.sam2=None
		self.out="out"
		self.zeroinfer=False
		self.ambigskip=False
		self.delim=None
		self.minlen=1

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "h" or opt == "help":
				continue
			elif opt=="1":
				self.sam1=arg
			elif opt=="2":
				self.sam2=arg
			elif opt=="r":
				self.ref=arg
			elif opt=="o":
				self.out=arg
			elif opt=="z":
				self.zeroinfer=True
			elif opt=="n":
				self.ambigskip=True
			elif opt=="d":
				self.delim=arg
			elif opt=="m":
				self.minlen=int(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.one2many and not self.many2one:
			self.display_help("No files provided.")



	def display_help(self, message=None):
		"""Display help menu"""
		if message is not None:
			print()
			print (message)
		print ("\nfastaFormatter.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Computes allele depth ratios from pileup data")
		print("""
	Mandatory arguments:
		-r	: Reference FASTA file
		-1	: Sample 1 pileup file
		-2	: Sample 2 pileup file
	Optional arguments:
		-z	: Infer zeros for positions missing in -1 and -2
			  NOTE: By default missing positions will be excluded from calculations
		-n	: Only count non-ambiguous (N) positions in reference
		-d	: FASTA header delimiter, after which characters are skipped [default=None]
		-m	: Minimum scaffold length to report [default=None]
		-o,--out	: Output file prefix [default=out]
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()

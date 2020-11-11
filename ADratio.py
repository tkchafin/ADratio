#!/usr/bin/python

import sys
import os
import getopt
import numpy as np
import pandas as pd
import naiveBayes as nb

def main():
	params = parseArgs()
	
	#build dictionary of reference contigs
	reference_lengths=dict()
	reference_ambigs=dict()
	maxn=0
	minlen=0
	read=0
	kept=0
	if params.ref:
		print("\nReading reference genome",params.ref)
		for contig in read_fasta(params.ref):
			read+=1
			if len(contig[1]) >= params.minlen:
				name=contig[0]
				if params.delim:
					name=(contig[0].split(params.delim))[0]

				#If -n, gather map of N positions
				pos=getSequencePositions(contig[1], "N", case_sensitive=False)
				
				#if too many Ns, skip contig
				#print("Number Ns:",len(pos))
				#print(float(len(pos))/float(len(contig[1])))
				if (float(len(pos))/float(len(contig[1]))) > params.maxambig:
					maxn+=1
				else:
				#else (only keep data for passing scaffolds)
				#keep ambiguous positions in a dictionary
					if params.ambigskip:
						if name not in reference_ambigs.keys():
							reference_ambigs[name] = list()
						reference_ambigs[name].append(pos)
					#also save reference lengths for fast lookup later
					reference_lengths[name] = len(contig[1])
					kept+=1
			else:
				minlen+=1
				continue
		print("\n\nTotal contigs read:",str(read))
		if minlen > 0:
			print("Contigs skipped below min length:",str(minlen))
		if maxn > 0:
			print("Contigs skipped above max N proportion:",str(maxn))
		print("Kept",str(kept),"contigs.\n")
	else:
		print("ERROR: No reference FASTA provided.")
		sys.exit()
	
	#Parse individual 1 coverage
	print("\nParsing bedgraph for individual 1:",params.sam1)
	if params.ambigskip:
		bad=reference_ambigs
	else:
		bad=None
	cov1 = parseBedGraph(params.sam1, reference_lengths, bad)
	print(cov1)
	
	#parse individual 2 coverage
	print("\nParsing bedgraph for individual 2:",params.sam2)
	if params.ambigskip:
		bad=reference_ambigs
	else:
		bad=None
	cov2 = parseBedGraph(params.sam2, reference_lengths, bad)
	print(cov2)
	
	#Calculate normalized ADratio per-scaffold
	print("\nComputing ADratio for each scaffold...")
	print("...Using the normalizing constant:",str(params.constant))
	adratios = computeADratios(cov1, cov2, params.constant)
	print(adratios)
	
	#if classification requested
	if params.classify:
		classifier = nb.nbClassifier()
		if params.config:
			priors = pd.read_table(params.config, sep="\t", header=0)
		else:
			priors = getDefaultPriors()
		classifier.initFromTable(priors)
		print("\nClassifying scaffolds using the following priors:")
		classifier.printPriors()
		
		dat=pd.read_table("example/nb_testdata.txt", sep="\t", header=0)
		dat=classifier.classify(dat)
		dat=classifier.getMAP(dat, params.map_thresh)
		if params.jaynes:
			dat=classifier.getJaynes(dat, params.j_thresh)
		with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
			print(dat)
		sys.exit()

def getDefaultPriors():
	d = {'Class': ["X", "Y", "auto"], 
	'AD_mean': [2.0, 0.0, 1.0],
	'AD_sd' : [0.1, 0.1, 0.1], 
	'Prob' : [1.0, 1.0, 1.0]}
	return(pd.DataFrame(data=d))

def computeADratios(cov1, cov2, constant=1.0):
	ad=dict()
	all = set(list(cov1.keys()) + list(cov2.keys()))
	for contig in all:
		if contig in cov1.keys():
			d1=float(cov1[contig])
		else:
			d1=0.0
		if contig in cov2.keys():
			d2=float(cov2[contig])
		else:
			d2=0.0
		ad[contig] = (d1/d2) * constant
	return(ad)

def parseBedGraph(bgfile, lens, bad_nucs=None):
	"""
	Function parses a bedgraph file contig-by-contig to calculate
	mean per-scaffold depths
	"""
	myCoverage = dict()
	this_scaffold=None
	if os.path.exists(bgfile):
		with open(bgfile, 'r') as fh:
			try:
				contig = ""
				seq = ""
				depths=None
				for line in fh:
					line = line.strip()
					if not line:
						continue
					stuff = line.split()
					if stuff[0] not in lens.keys():
						continue
					elif this_scaffold and this_scaffold != stuff[0]:
						#calculate depth of previous scaffold and set this_scaffold to current one
						if bad_nucs and this_scaffold in bad_nucs.keys():
								bad=bad_nucs[this_scaffold]
						else:
							bad=None
						myCoverage[this_scaffold] = getMeanDepth(depths, bad)
						depths=np.zeros(lens[this_scaffold])
						this_scaffold=stuff[0]
					elif not this_scaffold:
						#must be a new scaffold
						this_scaffold=stuff[0]
						if this_scaffold not in lens.keys():
							continue
						depths = np.zeros(lens[this_scaffold])
					start=int(stuff[1])
					end=int(stuff[2])
					cov=int(stuff[3])
					if(end-1 >= lens[this_scaffold]):
						print("Warning: Interval exceeds length of contig:",line)
					depths[start:end] = cov
				#add final scaffold and return dict
				myCoverage[this_scaffold] = getMeanDepth(depths, bad)
				return(myCoverage)
			except IOError:
				print("Could not read file ",bgfile)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%bgfile)	
		
#get mean from np.array of depths
def getMeanDepth(depths, exclude_chars=None):
	"""
	Function computes the mean from a np.array of integers
	        representing per-base sequencing depths. 
	Inputs: np.array of integers, (Optional) list of array indices
	        to exclude from computations
	Output: mean of input array, excluding nan values and exclude_chars
	"""
	#print(depths)
	if exclude_chars:
		depths[tuple(exclude_chars)] = np.nan
	#print(depths)
	mean_depth = np.nanmean(depths)
	return(mean_depth)

#Return indexes in string matching requested character
def getSequencePositions(seq, char, case_sensitive=True):
	"""
	Function returns all indices of a given character in a string, 
	allowing for either case sensitive or case insensitive comparisons
	Inputs: A string, a character, and (optional) case_sensitive boolean
	Outputs: List of (integer) indices 
	"""
	if not case_sensitive:
		#print(seq)
		char=char.lower()
		return([pos for pos, c in enumerate(seq) if c.lower() == char])
	else:
		return([pos for pos, c in enumerate(seq) if c == char])

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
			options, remainder = getopt.getopt(sys.argv[1:], 'h1:2:r:o:znd:m:M:c:bNp:P:xJj:', \
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
		#self.zeroinfer=False
		self.ambigskip=False
		self.delim=None
		self.minlen=1
		self.bedgraph=True #NOTE: Might add support for other formats later
		self.maxambig=0.5
		self.constant=1.0
		self.classify=False
		self.config=None
		self.map_thresh=0.0
		self.noPlots=False
		self.jaynes=False
		self.j_thresh=30

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
			elif opt=="b":
				self.bedgraph=True
			elif opt=="c":
				self.constant=float(arg)
			elif opt=="M":
				self.maxambig=float(arg)
			elif opt=="x":
				self.noPlots=True
			elif opt=="N":
				self.classify=True
			elif opt=="p":
				self.config=arg
			elif opt=="P":
				self.map_thresh=float(arg)
			elif opt=="J":
				self.jaynes=True
			elif opt=="j":
				self.j_thresh=float(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.ref:
			self.display_help("No reference FASTA provided.")
		if not self.sam1 or not self.sam2:
			self.display_help("Must provide 2 sample files (-1, -2)")



	def display_help(self, message=None):
		"""Display help menu"""
		if message is not None:
			print()
			print (message)
		print ("\nADratio.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Computes allele depth ratios from pileup data")
		print("""
	Mandatory arguments:
		-r	: Reference FASTA file
		-1	: Sample 1 coverage file
		-2	: Sample 2 coverage file
		
	Optional arguments:
		-d	: FASTA header delimiter [default=None]
		-o	: Output file prefix [default=out]
		
	ADratio arguments:
		-c	: Normalizing constant, calculated as:
			  # Sample 1 reads / # Sample 2 reads [Default=1.0]
		-n	: Only count non-ambiguous (N) positions in reference
		-m	: Minimum scaffold length to report [default=None]
		-M	: Maximum proportion of Ns to retain a contig [default=0.5]
		
	Classifier arguments:
		-N	: Classify scaffolds to chromosome type (e.g. X, Y, autosome)
		-p	: (Optional) Params file to customize chr type priors
			   See documentation. By default, we assume three Gaussian 
			   priors representing how we expect ADratio to vary by chr type:
			   Class	AD_mean	AD_sd	Prob
			   X	2.0	0.1	1.0
			   Y	0.0	0.1	1.0
			   auto	1.0	0.1	1.0
			   NOTE: Here we assume <-1> female and <-2> male.
		-P	: Maximum a posteriori threshold to keep a classification 
		-J	: Toggle on to calculate Jayne's 'evidence' for each class
		-j	: Jayne's evidence (db) threshold [default=30]
		
	Plotting arguments:
		-x	: Toggle to turn OFF plotting
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()

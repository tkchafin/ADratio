#!/usr/bin/python

import sys
import os
import getopt
import numpy as np
import pandas as pd
import nbClassifier as nb
import plotAD as plot

def main():
	params = parseArgs()
	
	#build dictionary of reference contigs
	reference_lengths=dict()
	reference_ambigs=dict()
	maxn=0
	minlen=0
	read=0
	kept=0
	if params.resume == 0:
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
			if kept <1:
				print("\nNothing to do.\n")
		else:
			print("ERROR: No reference FASTA provided.")
			sys.exit()
	else:
		print("Resuming from existing files...\n")
	
	if params.resume != 1:
		#Parse individual 1 coverage
		print("\nParsing bedgraph for individual 1:",params.sam1)
		if params.ambigskip:
			bad=reference_ambigs
		else:
			bad=None
		cov1 = parseBedGraph(params.sam1, reference_lengths, bad)
		oo=params.out+"_ind1"
		writeCov(oo, cov1)
		#print(cov1)
		
		#parse individual 2 coverage
		print("\nParsing bedgraph for individual 2:",params.sam2)
		if params.ambigskip:
			bad=reference_ambigs
		else:
			bad=None
		cov2 = parseBedGraph(params.sam2, reference_lengths, bad)
		op=params.out+"_ind2"
		writeCov(op, cov2)
		#print(cov2)
	else:
		print("Reading coverage files...\n")
		o1=params.out + "_ind1_cov.txt"
		cov1=readCoverage(o1)
		o2=params.out + "_ind2_cov.txt"
		cov2=readCoverage(o2)
	
	if params.resume not in [2, 3]:
		#Calculate normalized ADratio per-scaffold
		print("\nComputing ADratio for each scaffold...")
		print("...Using the normalizing constant:",str(params.constant))
		adratios = computeADratios(cov1, cov2, params.constant)
		#print(adratios)
		dat=ADtoDF(adratios)
		writeAD(adratios, params.out)
		#delete cov data; not needed any more 
		del cov1
		del cov2
		del adratios
	else:
		print("Reading ADratio file...\n")
		oa=params.out + "_AD.txt"
		adratios=pd.read_table(oa, sep="\t", header=0)
		dat=ADtoDF(adratios)
		del adratios
	
	
	#if classification requested
	if params.resume != 3:
		if params.classify:
			classifier = nb.nbClassifier()
			if params.config:
				priors = pd.read_table(params.config, sep="\t", header=0)
			else:
				priors = getDefaultPriors()
			if not params.fitdata:
				classifier.initFromTable(priors)
			else:
				print("\nFitting classifier using provided dataset:",params.fitdata)
				fdf=pd.read_table(params.fitdata, sep="\t", header=0)
				classifier.fit(fdf)
				if params.equalprobs:
					classifier.setProbsEqual()
			print("\nClassifying scaffolds using the following priors:")
			classifier.printPriors()
			
			#dat=pd.read_table("example/nb_testdata.txt", sep="\t", header=0)
			#dat=pd.read_table("example/nb_testfit.txt", sep="\t")
			dat_class=classifier.classify(dat)
			dat_class=classifier.getMAP(dat_class, params.map_thresh)
			if params.jaynes:
				dat_class=classifier.getJaynes(dat_class, params.j_thresh)
			#prettyPrint(dat_class)
			
			#output classification table
			oname=params.out+"_classify.txt"
			print("\nOutputting classified results to:",oname)
			dat_class.to_csv(oname, sep="\t", header=True, quoting=None, index=False)
		
	#Make plots
	if not params.noPlots:
		#histogram of AD values
		if params.classify:
			o=params.out + "_MAP"
			plot.plotADclassified(dat_class, o, "MAP", stat=params.stat, binwidth=params.binwidth, y=params.ylim, x=params.xlim)
			if params.jaynes:
				o2=params.out + "_JAYNE"
				plot.plotADclassified(dat_class, o2, "JAYNE", stat=params.stat, binwidth=params.binwidth, y=params.ylim, x=params.xlim)
		else:
			#dat=pd.read_table("example/nb_testfit.txt", sep="\t")
			plot.plotAD(dat, params.out, stat=params.stat, binwidth=params.binwidth, y=params.ylim, x=params.xlim)
			
		
	print("\nDone!\n")
	sys.exit()

def writeCov(out, cov):
	oname=out+"_cov.txt"
	new=dict()
	new["Scaffold"] = list(cov.keys())
	new["MeanDepth"] = list(cov.values())
	print("Outputting mean coverages to:",oname)
	df=pd.DataFrame.from_dict(new, orient="columns")
	df.to_csv(oname, sep="\t", header=True, quoting=None, index=False)

def readCoverage(f):
	df=pd.read_table(f, sep="\t", header=0)
	return(dict(zip(df.Scaffold,df.MeanDepth)))

def writeAD(ad, out):
	oname=out + "_AD.txt"
	df=ADtoDF(ad)
	print("Outputting AD results to:",oname)
	df.to_csv(oname, sep="\t", header=True, quoting=None, index=False)

def ADtoDF(ad):
	new=dict()
	new["Scaffold"]=list(ad.keys())
	new["AD"]=list(ad.values())
	df=pd.DataFrame.from_dict(new, orient="columns")
	#prettyPrint(df)
	return(df)

def prettyPrint(df):
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
		print(df.sort_values(df.columns[0], ascending=True))

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
			continue
		if d2 <= 0.0:
			continue
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
		#print(exclude_chars[0])
		j=sorted(exclude_chars[0], reverse=True)
		j2=[i for i in j if i < len(depths)]
		depths[j2] = np.nan
	#print(depths)
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
			options, remainder = getopt.getopt(sys.argv[1:], 'h1:2:r:o:znd:m:M:c:b:Np:P:xJj:F:fR:X:Y:S:', \
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
		self.binwidth=0.1
		self.maxambig=0.5
		self.constant=1.0
		self.classify=False
		self.config=None
		self.map_thresh=0.0
		self.noPlots=False
		self.jaynes=False
		self.j_thresh=30
		self.fitdata=None
		self.equalprobs=False
		
		self.resume=0
		self.xlim=None
		self.ylim=None
		self.stat='frequency'

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
				self.binwidth=float(arg)
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
			elif opt=="F":
				self.fitdata=arg
			elif opt=="f":
				self.equalprobs=True
			elif opt=="X":
				self.xlim=float(arg)
			elif opt=="Y":
				self.ylim=float(arg)
			elif opt=="S":
				if arg != "frequency" and arg != "count" and arg != "density" and arg != "probability":
					self.display_help("Invalid option for -S:",arg)
				else:
					self.stat=arg
			elif opt=="R":
				if int(arg) not in [1, 2, 3]:
					self.display_help("Invalid option for -R:",arg)
				else:
					self.resume=int(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.ref:
			self.display_help("No reference FASTA provided.")
		if not self.sam1 or not self.sam2:
			if self.resume == 0:
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
		-R	: Resume from: 1-Coverage files; 2-ADratio file; 3-Classify file
		
	ADratio arguments:
		-c	: Normalizing constant, calculated as:
			  # Sample 1 reads / # Sample 2 reads [Default=1.0]
		-n	: Only count non-ambiguous (N) positions in reference
		-m	: Minimum scaffold length to report [default=None]
		-M	: Maximum proportion of Ns to retain a contig [default=0.5]
		
	Classifier arguments:
		-N	: Toggle to classify scaffolds to chromosome type (e.g. X, Y, aut)
		-p	: (Optional) Params file to customize chr type priors
			   See documentation. By default, we assume three Gaussian 
			   priors representing how we expect ADratio to vary by chr type:
			   Class	AD_mean	AD_sd	Prob
			   X	2.0	0.1	1.0
			   Y	0.0	0.1	1.0
			   auto	1.0	0.1	1.0
			   NOTE: Here we assume <-1> female and <-2> male.
		-F	: (Optional) Fit classifier to a tab-delimited data file
			   See documentation. Format should be like so: 
			   Class AD
			   X	2.09
			   X	1.99
			   Y	0.001
			   auto	1.1
			   auto	0.977
			   ...
		-f	: Toggle on to set class probabilities equal when using -F
		-P	: Maximum a posteriori (MAP) threshold to keep a classification 
		-J	: Toggle on to calculate Jayne's 'evidence' for each class
		-j	: Jayne's evidence (db) threshold [default=30]
		
	Plotting arguments
		-x	: Toggle to turn OFF plotting
		-b	: Binwidth for plotting [default=0.1]
		-X	: X-limit for plotting [default=None]
		-Y	: Y-limit for plotting [default=None]
		-S	: Histogram stat [default='frequency']
			  Options: count, frequency, density, probability
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()

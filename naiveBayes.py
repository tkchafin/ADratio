#!/usr/bin/python

"""
Special thanks to tutorial by Jason Brownlee at Machine Learning Mastery:
https://machinelearningmastery.com/naive-bayes-classifier-scratch-python/
"""

import sys
import os
import re
import getopt
import numpy as np
import pandas as pd
#from scipy.stats import norm
import collections
from math import sqrt
from math import exp
from math import pi
from collections import OrderedDict

#nbClassifier Prior; consists of a Gaussian mean and sd
class nbPrior():
	"""
	Object holding Gaussian prior for a single variable
	Member of nbClass.priors
	"""
	def __init__(self, mean, sd):
		self.mean=mean
		self.sd=sd
	
	def getMean(self):
		return(self.mean)
	
	def getSD(self):
		return(self.sd)
	
	def getProb(self, value):
		exponent = exp(-((value-self.mean)**2 / (2 * self.sd**2 )))
		return ((1 / (sqrt(2 * pi) * self.sd)) * exponent)
		#d=norm(loc=self.mean, scale=self.sd)
		#return(d.pdf(value))

#nbClassifier Class; each holds a list of priors
class nbClass():
	"""
	Member of nbClassifier.class. 
	Contains a set of Gaussian priors for a named class in nbClassifier
	"""
	def __init__(self, name, classProb, priors=None):
		if not priors:
			self.name=name
			if priors:
				self.priors=priors
			else:
				self.priors=OrderedDict()
		self.prob=classProb

	def addPrior(self, priorName, mean, sd):
		self.priors[priorName] = nbPrior(mean, sd)
	
	def getName(self):
		return(self.name)
	
	def getClassProb(self):
		return(self.prob)
	
	def setClassProb(self, value):
		self.prob=value
	
	def calculatePosterior(self, values):
		prob = self.prob
		for index, v in enumerate(values):
			#multiply prob of each variable given prior
			prob *= list(self.priors.items())[index][1].getProb(v)
			return(prob)
	

#simple naive Bayes classifier
class nbClassifier():
	"""
	Class provides text command-line menu and holds parameters
	needed for the run
	"""
	def __init__(self):
		self.classes=OrderedDict()
		self.nClasses=0
		self.priors=set() #keep track of prior names
		self.nPriors=0
		self.classProbs=list()
	
	def addClass(self, className, classProb, priors=None):
		"""
		Adds a class to the classifier, with or without specified priors
		"""
		self.classes[className] = nbClass(className, classProb, priors)
		self.classProbs.append(classProb)
		self.nClasses += 1
	
	def addPrior(self, className, priorName, mean, std):
		"""
		Adds a Gaussian prior to a class contained in the classifier
		"""
		if className in self.classes:
			self.classes[className].addPrior(priorName, mean, std)
			if priorName not in self.priors:
				self.priors.add(priorName)
				self.nPriors += 1
		else:
			raise Exception("nbClassifier, no class by name",className)
	
	def initFromTable(self, table):
		"""
		Function manually specifies priors for each class from a 
		pandas dataframe of the form:
		class   varA_mean    varA_sd    varB_mean .....
		"""
		#print(table)
		variables=[ x.split("_")[0] for x in list(table.columns[:-1])[1::2]]
		#print(variables)
		for index, row in table.iterrows():
			#print(row)
			self.addClass(str(row[0]), float(row[-1]))
			#self.classProbs.append(float(row[-1]))
			offset=0
			#print(row)
			for v in variables:
				self.addPrior(str(row[0]), v, float(row[offset+1]), float(row[offset+2]))
				offset += 2
	
	def printPriors(self):
		"""
		Prints a table of Gaussian priors and P(class) for each class.
		"""
		d=dict()
		d['Class']=list(self.classes.keys())
		for v in self.priors:
			m=str(v) + "_mean"
			s=str(v) + "_sd"
			d[m] = list()
			d[s] = list()
			for c in d['Class']:
				d[m].append(self.classes[c].priors[v].getMean())
				d[s].append(self.classes[c].priors[v].getSD())
		d['classProb'] = self.classProbs
		#print(d)
		df = pd.DataFrame(data=d)
		with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
			print(df)
	
	def classify(self, data):
		"""
		Function applies a fitted naive bayes classifier to 
		a set of observation data. Observations should be formatted
		as a pandas dataframe, where the first column is the record name,
		and following columns correspond to the class prior variable names. 
		Order does not matter. 
		
		Function returns a modified dataframe with a columns for each CLASS,
		and values corresponding to probabilities of assignment
		"""
		#get indices of each variable in dataset
		indices = list()
		for v in self.priors:
			indices.append(list(data.columns).index(v))
		
		#for each row, compute probability of membership for each class
		probs=OrderedDict()
		for index, row in data.iterrows():
			name=row[0]
			#print("Name:",name)
			values=[row[i] for i in indices]
			#print("Data:",values)
			for k in self.classes:
				c=self.classes[k].getName()
				if c not in probs.keys():
					probs[c] = list()
				p=self.classes[k].calculatePosterior(values)
				#print(p)
				probs[c].append(p)
		#add probabilities to table
		for k in probs.keys():
			data[k] = probs[k]
		return(data)	
	
	def getMAP(self, data, threshold=0.0):
		"""
		Given posterior probabilities for each class, 
		function adds two columns to dataframe:
		1) MAP_value: The maximum a posteriori estimate
		2) Class corresponding to that estimate.
		
		Also accepts threshold>0.0 argument, which limits classifications to 
		those above a probability threshold
		"""
		MAP_values = data[list(self.classes.keys())].max(axis=1)
		MAP_names = data[list(self.classes.keys())].idxmax(axis=1)
		data["MAP_value"] = MAP_values
		data["MAP"] = MAP_names
		if threshold > 0.0:
			data.loc[data['MAP_value'] < threshold, ['MAP']] = "None"
		return(data)
	
	def getJaynes(self, data, threshold=0.0):
		"""
		Function computes Jayne's evidence from class assignment
		posterior probabilities as J = prob(class) / product(prob(other_classes)),
		then transformed as 10*log10(J) for interpretability
		
		Adds two columns to dataframe:
		1) JAYNE_value: The highest evidence value among all classes
		2) JAYNE: The class with the highest evidence
		
		Also accepts threshold>0.0 argument, which limits classifications to 
		those above a evidence threshold
		"""
		#get indices for each class prob
		indices = list()
		for v in self.priors:
			indices.append(list(data.columns).index(v))
		#for each class, get prob(class) / prob(not-class)
		j = dict()
		j_vars=list()
		for k in self.classes.keys():
			p_class = data[k]
			p_notclass = 1.0
			for kk in self.classes.keys():
				if kk != k:
					p_notclass *= data[kk]
			p_ratio = p_class/p_notclass
			j_class = 10 * np.log10(p_ratio)
			name = str(k) + "_J"
			data[name] = j_class
			j_vars.append(name)
			
		J_values = data[j_vars].max(axis=1)
		J_names = data[j_vars].idxmax(axis=1)
		data["JAYNE_value"] = J_values
		data["JAYNE"] = J_names
		data["JAYNE"] =  [re.sub(r'_J','', str(x)) for x in data["JAYNE"]]
		
		#remove classification if below evidence threshold
		if threshold > 0.0:
			data.loc[data['JAYNE_value'] < threshold, ['JAYNE']] = "None"
		
		return(data)
		
	def fit(self, data):	
		"""
		Function fits a Naive Bayes classifier to a given dataset
		Dataset should consist of N columns where column 1 
		is the "class" (or group) assignment of the entry, 
		and all columns 2-N are values. 
		"""
		variables=list(data.columns)[1:]
		total=data.shape[0]
		for name, group in data.groupby(data.columns[0]):
			class_num = group.shape[0]
			class_prob= class_num/total
			self.addClass(name, class_prob)
			for var in variables:
				var_mean = np.nanmean(group[var])
				var_sd = np.nanstd(group[var])
				#print(name, var, var_mean, var_sd, class_prob)
				self.addPrior(name, var, var_mean, var_sd)
	
	def setProbsEqual(self):
		"""
		function sets all class probabilities equal (to 1.0)
		"""
		idx=0
		for c in self.classes.keys():
			self.classes[c].setClassProb(1.0)
			self.classProbs[idx]=1.0
			idx+=1
	
	def setClassProbs(self, l):
		"""
		set class probabilitis given a vector of values
		"""
		idx=0
		for c in self.classes.keys():
			self.classes[c].setClassProb(float(l[idx]))
			self.classProbs[idx]=1.0
			idx+=1
					
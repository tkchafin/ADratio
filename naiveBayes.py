#!/usr/bin/python

"""
Special thanks to tutorial by Jason Brownlee at Machine Learning Mastery,
from which some of this code was taken. 
https://machinelearningmastery.com/naive-bayes-classifier-scratch-python/
"""

import sys
import os
import re
import getopt
import numpy as np
import pandas as pd
import collections
from math import sqrt
from math import exp
from math import pi
from collections import OrderedDict

#nbClassifier Prior; consists of a Gaussian mean and sd
class nbPrior():
	def __init__(self, mean, sd):
		self.mean=mean
		self.sd=sd
	
	def getMean(self):
		return(self.mean)
	
	def getSD(self):
		return(self.sd)
	
	def getProb(self, value):
		exponent = exp(-((value-self.mean)**2 / (2 * self.sd**2 )))
		return (1 / (sqrt(2 * pi) * self.sd)) * exponent

#nbClassifier Class; each holds a list of priors
class nbClass():
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
		return(self.classProb)
	
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
		self.classes[className] = nbClass(className, classProb,priors)
		self.nClasses += 1
	
	def addPrior(self, className, priorName, mean, std):
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
			self.classProbs.append(float(row[-1]))
			offset=0
			print(row)
			for v in variables:
				self.addPrior(str(row[0]), v, float(row[offset+1]), float(row[offset+2]))
				offset += 2
	
	def printPriors(self):
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
		df = pd.DataFrame(data=d)
		with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
			print(df)
	
	def classify(self, data):
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
	
	def getMAP(self, data, threshold):
		MAP_values = data[list(self.classes.keys())].max(axis=1)
		MAP_names = data[list(self.classes.keys())].idxmax(axis=1)
		data["MAP_value"] = MAP_values
		data["MAP"] = MAP_names
		if threshold > 0.0:
			data.loc[data['MAP_value'] < threshold, ['MAP']] = "None"
		return(data)
	
	def getJaynes(self, data, threshold=0.0):
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
		pass


# Calculate the mean of a list of numbers
def mean(l):
	return sum(l)/float(len(l))

# Calculate the standard deviation of a list of numbers
def stdev(l):
	avg = mean(l)
	var = sum([(x-avg)**2 for x in l]) / float(len(l)-1)
	return sqrt(var)		
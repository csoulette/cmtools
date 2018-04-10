#!/usr/bin/env python3


########################################################################
# File: 
#  executable: 
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 03/22/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
import re

from multiprocessing import Pool
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
from samJuncs import SAM
import pysam

########################################################################
# CommandLine
########################################################################

class CommandLine(object) :
	'''
	Handle the command line, usage and help requests.
	CommandLine uses argparse, now standard in 2.7 and beyond. 
	it implements a standard command line argument parser with various argument options,
	and a standard usage and help,
	attributes:
	myCommandLine.args is a dictionary which includes each of the available command line arguments as
	myCommandLine.args['option'] 
	
	methods:
	
	'''
	
	def __init__(self, inOpts=None) :
		'''
		CommandLine constructor.
		Implements a parser to interpret the command line argv string using argparse.
		'''
		import argparse
		self.parser = argparse.ArgumentParser(description = ' TBD - TBD.',
											 epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
											 add_help = True, #default is True 
											 prefix_chars = '-', 
											 usage = '%(prog)s -i aligned.bam ')
		# Add args
		self.parser.add_argument('-i', "--input_bam", action = 'store', required=True, help='Input bam file.')
		self.parser.add_argument('-j1', "--annotated_junctions", action = 'store', required=False, help='Input bed file.')
		self.parser.add_argument('-j2', "--other_junctions", action = 'store', required=False, help='Input bed file.')
		self.parser.add_argument('-tss', "--known_starts", action = 'store', required=False, help='Input bed file.')
		self.parser.add_argument('-tse', "--known_ends", action = 'store', required=False, help='Input bed file.')
		self.parser.add_argument('--isSam', action = 'store_true',  default=False, help='File is sam')
		self.parser.add_argument('-p', "--threads", action = 'store', required=False, default = 2, help='Input bed file.')
		self.parser.add_argument("--fix_jcn", action = 'store', required=False, default = False, help='Fix junctions using the wiggle rule (in development).')

		if inOpts is None :
			self.args = vars(self.parser.parse_args())
		else :
			self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Functions
########################################################################

def resolveJcns(chrom, j1, j2, real, jcns, jcnTree):

	try:
		leftHits = [x.data for x in jcnTree[chrom][j1]]

		bestLeft = min(leftHits, key=lambda x:abs(x-j1))

		rightHits = [x.data for x in jcnTree[chrom][j2]]

		bestRight = min(rightHits, key=lambda x:abs(x-j2))

		return(bestLeft,bestRight)

	except:
		#No best match
		return (None, None)

def jcnDB(known, other, wiggle):

	# First define set of junctions to be quantified
	realJcn = dict()
	realSites = dict()
	siteTrees = dict()
	wig = 10
	with open(known, 'r') as lines:
		for line in lines:
			chrom, c1, c2, txn, score, strand = line.rstrip().split()

			if chrom not in realJcn:
				realJcn[chrom] = dict()
				siteTrees[chrom] = IntervalTree()
				realSites[chrom] = dict()

			c1, c2 = int(c1), int(c2)

			jcn = (c1, c2)
			realJcn[(c1,c2)] = "K"
			realSites[c1] = "K"
			realSites[c2] = "K"
			siteTrees[chrom][c1-wig:c1+wig] = c1
			siteTrees[chrom][c2-wig:c2+wig] = c2

	with open(other, 'r') as lines:
		for line in lines:
			chrom, c1, c2, scores = line.rstrip().split()
			c1, c2 = int(c1), int(c2)

			if c1 in realSites and c2 in realSites:
				continue

			elif c1 not in realSites and c2 not in realSites:
				realSites[c2] = "S"
				realSites[c1] = "S"
				realJcn[(c1,c2)] = "S"
				siteTrees[chrom][c1-wig:c1+wig] = c1
				siteTrees[chrom][c2-wig:c2+wig] = c2

			elif c1 not in realSites:
				strand = realSites[c2]
				realSites[c1] = "S"
				realJcn[(c1,c2)] = "S"
				siteTrees[chrom][c1-wig:c1+wig] = c1

			elif c2 not in realSites:
				strand = realSites[c1]
				realSites[c2] = "S"
				realJcn[(c1,c2)] = "S"
				siteTrees[chrom][c2-wig:c2+wig] = c2

	return [realJcn, realSites, siteTrees]

########################################################################
# Main
# Here is the main program
# 
########################################################################
def main():
	'''
	stuff...
	'''

	myCommandLine = CommandLine()
	
	alignmentFile = myCommandLine.args['input_bam']
	knownSites = myCommandLine.args['annotated_junctions']
	otherSites = myCommandLine.args['other_junctions']
	starts = myCommandLine.args['known_starts']
	ends = myCommandLine.args['known_ends']
	threads = myCommandLine.args['threads']
	isBam = myCommandLine.args['isSam']	
	fix = myCommandLine.args['fix_jcn']

	# First define set of junctions to be quantified
	#realJcn, realSites, siteTrees = jcnDB(knownSites, otherSites, 10)

	sObj = SAM(alignmentFile, isBam)

	# TO do multithreading
	#headers = re.findall('SN:(\S+)', pysam.view(alignmentFile, "-H"))
	#print(headers)
	#sys.exit(1)
	#p = Pool(threads)

	if fix:
		for num, readData in enumerate(sObj.readJuncs(),0):
			#do something
			read, chrom, startPos, junctions, endPos = readData
		
			#for jcn in list(junctions:

			print(chrom, startPos, endPos, read, sep="\t")
	else:
		for num, readData in enumerate(sObj.readJuncs(),0):
			read, chrom, startPos, junctions, endPos, flags, tags = readData

			print(chrom, startPos, endPos, read, ".", flags, "\t".join(["%s:%s" % (x[0],x[1]) for x in junctions]), tags, sep="\t")			

if __name__ == "__main__":
	main();        

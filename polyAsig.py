from __future__ import print_function

########################################################################
# File: polyAsig.py
#  executable: polyAsig.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 01/17/2019 Created
#
########################################################################

import sys, os
from Bio import motifs
import re
import random
random.seed(123)
fasta = sys.argv[1]
seqs = list()


signal1 = "AATAAA"[::-1]
signal2 = "ATTAAA"[::-1]
signal3 = "AGTAAA"[::-1]
random4 = "".join([random.choice('TCAG') for _ in range(6)])[::-1]

with open(fasta) as lines:
	for header in lines:
		sequence = next(lines).rstrip()
		if len(sequence)<150: continue
		endChunk = sequence[-150:]
		endChunk = endChunk[::-1]
		seqs.append(endChunk)
		try:
			print(re.search(signal1, endChunk).start(),signal1[::-1])
			print(re.search(signal2, endChunk).start(),signal2[::-1])
			print(re.search(signal3, endChunk).start(),signal3[::-1])
			print(re.search(random4, endChunk).start(),random4[::-1])
		except:
			
			continue

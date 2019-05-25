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


fasta = sys.argv[1]

with open(fasta) as lines:
	for header in lines:
		sequence = next(lines)
		if len(sequence)<150: continue
		endChunk = sequence[-100:]

		print(endChunk)
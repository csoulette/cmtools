from __future__ import print_function

########################################################################
# File: coa.py
#  executable: coa.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 05/24/2019 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
from collections import Counter
import itertools
import re

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
        self.parser = argparse.ArgumentParser(description = ' tbd',
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -help')
        # Add args
        
        self.parser.add_argument("-i",'--isoforms',   action = 'store', required=True, help='FLAIR Isoforms')
        self.parser.add_argument("-c",'--counts_matrix',   action = 'store', required=True, help='counts matrix')

       
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

# gene

class Gene(object):
    '''
    '''
    def __init__(self, name = None):
        self.name = name
        self.isoforms = set()
        self.exonFeatures = set()
    
class Isoform(object):

    def __init__(self, name=None):
        self.name = name
        self.quants = None
        self.features = None
        self.exonSet = None

class Feature(object):

    def __init__(self, ftype=None, name=None):
        self.featureType = ftype
        self.name = name
        self.parents = set()

# funktions

def bed12toExons(start,starts,sizes):
    '''
    Take bed12 entry and convert block/sizes to exon coordinates.
    '''
    start = int(start)
    sizes, starts = list(map(int,sizes)), list(map(int,starts))
    exons = list()
    for num, st in enumerate(starts,0):
        c1 = st + start
        c2 = c1 + sizes[num]
        exons.append((c1,c2))
    return exons

def buildMatrices(counts, flairIsoforms):
    '''
    '''
    d = dict()
    features = dict()
    with open(flairIsoforms) as lines:
        for i in lines:
            cols = i.rstrip().split()
            starts,sizes = cols[-1].rstrip(",").split(","),cols[-2].rstrip(",").split(",")
            exons = bed12toExons(cols[1],starts,sizes)
            gene = re.search("(ENSG[^\.]+)", i).group(1)
            iso = cols[3]

            tempFeatureList = [] 
            for num, feature in enumerate(exons,1):
                if num == 0:
                    pass
                elif num == len(features):
                    pass
                else:
                    if feature not in features:
                        features[feature] = Feature("exon",feature)
                    feature = features[feature]
                    feature.parents.add(counts[iso])
                    tempFeatureList.append(feature)

            if gene not in d:
                d[gene] = Gene(gene)

            d[gene].exonFeatures = d[gene].exonFeatures.union(set(tempFeatureList))
            d[gene].isoforms.add(counts[iso])

    return d, counts


def parseCountsMatrix(counts):
    '''
    '''

    data = dict()
    with open(counts) as f:
        header = next(f)
        for i in f:
            cols = i.rstrip().split()
            iso, vals = cols[0],np.asarray(cols[1:], dtype=float)
            isoObj = Isoform(iso)
            isoObj.quants = vals
            data[iso] = isoObj
    return data


def getScores(counts,genes):

    for g,gobj in genes.items():
        
        combos = itertools.combinations(gobj.exonFeatures,2)
        for x in combos:
            left, right = x
            both = left.parents.intersection(right.parents)
            leftOnly = left.parents.difference(right.parents)
            rightOnly = right.parents.difference(left.parents)
            neither = gobj.isoforms.difference(left.parents.union(right.parents))

            bothQ = sum(np.asarray([x.quants for x in both]))
            leftQ = sum(np.asarray([x.quants for x in leftOnly]))
            rightQ = sum(np.asarray([x.quants for x in rightOnly]))
            neitherQ = sum(np.asarray([x.quants for x in neither]))
            print(left.name,right.name, bothQ, leftQ, rightQ, neitherQ, sep="\t")
            sys.exit(1)


def main():
    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()
    bed           = myCommandLine.args["isoforms"]
    counts        = myCommandLine.args["counts_matrix"]



    counts   = parseCountsMatrix(counts)
    genes, counts    = buildMatrices(counts, bed)
    print(genes)
    getScores(counts,genes)

  
# maine
if __name__ == "__main__":
    main()
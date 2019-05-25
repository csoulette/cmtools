from __future__ import print_function


########################################################################
# File: predictProductivity.py
#  executable: predictProductivity.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 10/09/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
from tqdm import *
import re
import pybedtools
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
        self.parser = argparse.ArgumentParser(description = ' predictProductivity - a tool.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -i isoforms.bed -f isoforms.fa -g annotations.gtf')
        # Add args
        self.parser.add_argument('-i', "--input_isoforms", action = 'store', required=True, help=' Input collapsed isoforms in bed12 format.')
        self.parser.add_argument('-g', "--gtf", action = 'store', required=True, help='Gencode annotation file.')
        self.parser.add_argument('-f', "--genome_fasta", action = 'store', required=True, help='Fasta file containing transcript sequences.')
        self.parser.add_argument("--quiet", action = 'store_false', required=False, default = True, help='Do not display progress')

        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Isoform
########################################################################

class Isoform(object) :
    '''
    Object to handle isoform related data.
    
    attributes:
        
    methods:
    
    ''' 

    def __init__(self, tid=None, gid=None, seq=None):
        self.tid = tid
        self.gid = gid
        self.pro = "UNK"
        self.chrom = ""

        self.sequence  = seq
        self.bed12     = None
        self.exons     = set()
        self.starts    = set()
        self.orfs      = list()
        self.orfStart  = int()
        self.orfEnd    = int()
        self.exonSizes = list()

    def sortORFs(self):

        self.orfs.sort(key=lambda x: x[-1])
        #for i in self.orfs:
        #    print(i,self.gid,self.tid)

########################################################################
# MAIN
########################################################################

def bed12ToExonRanges(cols):
    pass

def getStarts(gtf):
    starts = list()
    with open(gtf) as lines:
        for l in lines:
            if l[0] == "#": continue
            cols = l.rstrip().split("\t")
            chrom, c1, c2, strand = cols[0], int(cols[3])-1, int(cols[4]), cols[6]
            if cols[2] == "start_codon":
                
                gene = re.search("(ENSG[^\"]+)", cols[-1]).group(1)
                
                starts.append((chrom,c1,c2,gene,".",strand))

           
    return starts


def getSeqs(bed, genome):


    isoDict = dict()
    bt = pybedtools.BedTool(bed)
    bt.sequence(fi=genome, tab=True, s=True, split=True, name=True)
    with open(bt.seqfn) as entries:
        for entry in entries:
            read,seq  = entry.split()
            
            data = read.split("_")
            iso,gene = data[0], data[-1]
            if iso not in isoDict:
                isoDict[read] = Isoform(iso,gene,seq)
    return isoDict


def getStartRelPos(genomicStartPos,exon, exons, isoObj):
    '''
    is handed a genomic position, the exon it occurs in, all exons,
    and returns the position relative to all exons
    '''
    exonNum = exons.index(exon)
    isoObj.exonSizes = [x[1]-x[0] for x in exons]

    # First get start position relative to transcript sequence.
    if isoObj.strand == "+":
        relativeStart = genomicStartPos - exons[exonNum][0] + sum([x for x in isoObj.exonSizes[:exonNum]])
    elif isoObj.strand == "-":
        relativeStart = len(isoObj.sequence) - (genomicStartPos - exons[exonNum][0] + sum([x for x in isoObj.exonSizes[:exonNum]])) - 3
        
    return relativeStart


def checkPTC(orfEndPos, exons, isoObj):
    '''
    takes a transcript sequence position, and list of exon sizes to detemine
    if that position occurs more than 55nucleotides away from a splice junction.
    ptc = True if yes, ptc = False if not.
    the genomic position is also reported.
    '''
    stopDistFromExon = None
    exonWithStop = None 
    ptc  = None
    genomicPos = int()
    distance   = 0

    if isoObj.strand  == "-": 
        isoObj.exonSizes = isoObj.exonSizes[::-1]
        exons = exons [::-1]

    for num,e in enumerate(isoObj.exonSizes,0):

        distance += e

        # if the stop codon is in the last exon, then not ptc.
        if num == len(isoObj.exonSizes)-1:
            ptc = False
            if exonWithStop == None:
                exonWithStop = num
                stopDistFromExon = distance - orfEndPos

        # if the distance is greater than the stop position, then check if
        # the difference in distance is more than 55nt
        # if yet then ptc = True
        # also, track which exon the stop codon is in to get genomic position
        elif orfEndPos<distance:
            distToJunc =  distance - orfEndPos
            if exonWithStop == None:
                exonWithStop = num
                stopDistFromExon = int(distToJunc)

            if distToJunc>55:
                ptc=True
                break


    exonsWithStop = exons[exonWithStop]
    left,right    = exonsWithStop

    if isoObj.exonSizes[exonWithStop] != right - left:
        print(isoObj.strand,exons,isoObj.exonSizes,exonWithStop,exonsWithStop,orfEndPos,ptc)

    

    genomicPos = right - stopDistFromExon if isoObj.strand == "+" else left + stopDistFromExon 

    return genomicPos, ptc


def predict(bed, starts, isoDict):

    bt = pybedtools.BedTool(bed)
    b6 = bt.bed6()
    st = pybedtools.BedTool(starts)

    bt_st = b6.intersect(st, s=True, split=True, wao=True)
    for intersection in bt_st:
        read   = intersection[3]
        #iso,gene = read.split("_")
        overlap  = intersection[-1]
        goStart  = int(intersection[-6])
        exonCoord = (int(intersection[1]),int(intersection[2]))
        isoDict[read].strand = intersection[5]
        isoDict[read].chrom = intersection[0]
        isoDict[read].exons.add(exonCoord)

        if overlap != "3":
            continue
        else:
            isoDict[read].starts.add((exonCoord,goStart))

    stops = set(['TAA','TGA','TAG'])
    for iso,o in isoDict.items():
        exons = list(o.exons)
        exons.sort()

        if len(o.starts)<1:
            o.orfs.append(["NGO", exons[0][0], exons[0][0], 0])

        else:
            for start in o.starts:
                exon,startPos = start
                relativeStart = getStartRelPos(startPos,exon,exons,o)
                fiveUTR,rest  = o.sequence[:relativeStart], o.sequence[relativeStart:]
                
                # Next find first stop codon
                for i in range(0, len(rest), 3):
                    codon = rest[i:i+3]
                    if rest[i:i+3] in stops:
                        break

                # i is the last position after going through all codons and breaking at a stop
                # is a stop was never reached then i should represent the last NT in the entire seq
                # therefore, i+3 should be longer than the entire potential orf is a stop was never reached.
                # lets call these nonstop, or nst for now.
                if i+3 >= len(rest):
                    o.orfs.append(["NST", startPos, exons[-1][-1] if o.strand == "+" else exons[0][0], 0])  

                #else if a stop was reached...
                else:
                    orfEndPos = len(fiveUTR)+i+3
                    distance = 0
                    genomicStopPos, ptc = checkPTC(orfEndPos, exons, o)
                    ptc = "PTC" if ptc else "PRO"
                    o.orfs.append([ptc, startPos, genomicStopPos, orfEndPos - relativeStart ])


        o.sortORFs()

    return isoDict

def main():
    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()
    bed    = myCommandLine.args['input_isoforms']
    genome = myCommandLine.args['genome_fasta']
    gtf    = myCommandLine.args['gtf']

    starts      = getStarts(gtf)
    isoformObjs = getSeqs(bed, genome)
    isoformObjs = predict(bed, starts, isoformObjs)

    beaut = {"PRO":"103,169,207", "PTC":"239,138,98", "NST":"0,0,0","NGO":"0,0,0"}

    with open(bed) as lines:
        for line in lines:
            bedCols = line.rstrip().split()
            isoObj = isoformObjs[bedCols[3]]
            pro,start,end,orfLen = isoObj.orfs[-1]
            bedCols[3] = "%s_%s" % (bedCols[3], pro)
            bedCols[8] = beaut[pro]

            if isoObj.strand == "+":
                bedCols[6],bedCols[7] = str(start),str(end)
            else:
                bedCols[7],bedCols[6] = str(start),str(end)
            print("\t".join(bedCols))


if __name__ == "__main__":
    main()

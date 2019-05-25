from __future__ import print_function

########################################################################
# File: diaFLAIR.py
#  executable: diaFLAIR.py
# Purpose: wrapper for Differential Isoform Analyses
#
#          
# Author: Cameron M. Soulette
# History:      cms 01/17/2019 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import pandas as pd
import numpy as np
from subprocess import Popen

scriptPath = os.path.realpath(__file__)
path = "/".join(scriptPath.split("/")[:-1])
runDE = path + "/" + "runDE.py"
runDU = path + "/" + "runDU.py"
runAS = path + "/" + "runAS.py"
runAP = path + "/" + "runAP.py"
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
        self.parser = argparse.ArgumentParser(description = ' deFLAIR.py - a rpy2 convenience tool to run DESeq2.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s --manifest manifest.txt --workingdir dir_name --outdir out_dir --filter N')
        # Add args
        self.parser.add_argument("--outDir"    , action = 'store', required=True, 
                                    help='Write to specified output directory.')
        self.parser.add_argument("--filter"    , action = 'store', required=False, default = 10, type=int,
                                    help='Isoforms with less than specified read count for either Condition A or B are filtered (Default: 10 reads)')
        self.parser.add_argument("--manifest"    , action = 'store', required=True, 
                                    help='Tab separated file containing count file path, condition, and batch labels.')
        self.parser.add_argument("--isoforms"    , action = 'store', required=True, 
                                    help='FLAIR Isoform BED12 file.')
                


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

    def __init__(self, tid=None, gid=None, pro=None):
        self.tid = tid
        self.gid = gid
        self.pro = pro

        self.isoExp  = list()
        self.geneExp = list()
        self.usage   = list()

        self.dgeAdjPval = float()
        self.dieAdjPval = float()
        self.diuAdjPval = float()

        self.dgeFC = float()
        self.dieFC = float()
        self.diuDU = float()

    def computeUsage(self):

        self.usage = [np.divide(iso,gene) for iso,gene in zip(isoExp,geneExp)]
        self.usage[np.isinf(self.usage)] = np.nan


def filesToDF(f, thresh, isoforms):
    

    
    data  = dict()
    isoQ  = dict()
    geneQ = dict()
    samples = [x[0].split("/")[-1] for x in f]
    for num,i in enumerate(f,0):
        fName,group,batch = i
        with open(fName) as l:
            for line in l:
                name, count = line.rstrip().split()
                count = int(count) 
                iso, gene   = name, name.split("_")[-1]
                
                if gene not in geneQ:
                    geneQ[gene] =  np.zeros(len(f))
                    isoforms[iso].geneExp = np.zeros(len(f))

                if iso not in isoQ:
                    isoQ[iso] = np.zeros(len(f))
                    isoforms[iso].isoExp = np.zeros(len(f))

                isoQ[iso][num]    = count
                geneQ[gene][num] += count

                isoforms[iso].isoExp[num]    = count
                isoforms[gene].geneExp[num] += count
                


    # Convert counts to numpy array, where rows are read counts.
    # Also, Indices is an array of isoform/gene IDs.
    isoIndices  = np.asarray(list(isoQ.keys()))
    isoValues   = np.asarray([isoQ[x] for x in isoIndices], dtype=int)
    geneIndices = np.asarray(list(geneQ.keys()))
    geneValues  = np.asarray([geneQ[x] for x in geneIndices], dtype=int)
    
    # Filter count array and indices rows by threshold. 
    isoFiltered         = isoValues[(np.min(isoValues[:,3:],axis=1) > thresh) |  (np.min(isoValues[:,:3],axis=1) > thresh)]
    isoFilteredIndices  = isoIndices[(np.min(isoValues[:,3:],axis=1) > thresh) |  (np.min(isoValues[:,:3],axis=1) > thresh)]
    geneFiltered        = geneValues[(np.min(geneValues[:,3:],axis=1) > thresh) |  (np.min(geneValues[:,:3],axis=1) > thresh)]
    geneFilteredIndices = geneIndices[(np.min(geneValues[:,3:],axis=1) > thresh) |  (np.min(geneValues[:,:3],axis=1) > thresh)]
    
    isoDF  = pd.DataFrame(isoFiltered,columns=samples)
    geneDF = pd.DataFrame(geneFiltered,columns=samples)
    
    isoDF['ids']  = isoFilteredIndices
    geneDF['ids'] = geneFilteredIndices
    
    isoDF  = isoDF.set_index('ids')
    geneDF = geneDF.set_index('ids')

    return isoDF, geneDF, samples, isoforms


def makeDir(out):
    try:
        os.mkdir("./%s" % out)
    except:
        #exists
        pass

def checkFile(f):

    if os.path.isfile(fname):
        return f
    else:
        print("Cannot find quantification file %s. Exiting." % f, file=sys.stderr)
        sys.exit(1)

def checkSamples(manifest):


    groups = set()
    batches = set()
    sampleData = list()

    with open(manifest) as lines:
        # group batch sampleID file
        for l in lines:
            col = l.rstrip().split()
            fname, group, batch = col
            if len(col)!=3:
                print("Manifest does not have 3 columns. Please check format. Exiting.", file=sys.stderr)
                sys.exit(1)
            groups.add(group)
            batches.add(batch)
            sampleData.append(col)

    if len(list(groups))!=2:
        print("Number of conditions/groups does not equal 2. Exiting", file=sys.stderr)
        sys.exit(1)
    return list(groups),list(batches),sampleData

def flairIsoformsToObject(flairBED):
    '''
    Takes bed12 file containing FLAIR isoforms and returns
    isoform dict with Isoform Objects.
    '''
    isos = dict()
    with open(flairBED) as line:
        for l in lines:
            cols = l.rstrip().split()
            try:
                tid,gid,pro = cols[3].split("_")
            except:
                #misformatted line
                continue
            isos[tid] = Isoform(tid,gid,pro)
    return isos


def main():

    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()

    outDir     = myCommandLine.args['outDir']
    manifest   = myCommandLine.args['manifest']
    sFilter    = myCommandLine.args['filter']
    flairBED   = myCommandLine.args['isoforms']


    # Create output directory.
    makeDir(outDir)

    # Generate isoform objects.
    isoforms = flairIsoformsToObject(flairBED)

    # Check manifest formatting.
    group, batches, sampleData = checkSamples(manifest)
    group,batches = list(group), list(batches)

    # Convert count tables to dataframe and update isoform objects.
    isoformDF, geneDF, samples, isoforms = filesToDF(sampleData, sFilter, isoforms)

    for i,o in isoforms.items():
        o.computeUsage()
        print(o.isoExp, o.usage)
    sys.exit(1)
    
    header        = ['sampleName','condition','batch']
    formulaMatrix = [[x[0].split("/")[-1],x[1],x[2]] for x in sampleData]
    formulaDF     = pd.DataFrame(formulaMatrix,columns=header)
    formulaDF     = formulaDF.set_index('sampleName')

    formulaMatrixFile = "./%s/formula_matrix.tsv" % outDir
    isoMatrixFile     = "./%s/isoform_quant_matrix_deFLAIR.tsv" % outDir
    geneMatrixFile    = "./%s/gene_quant_matrix_deFLAIR.tsv" % outDir
    
    

    formulaDF.to_csv( formulaMatrixFile, sep='\t')
    isoformDF.to_csv( isoMatrixFile, sep='\t')
    geneDF.to_csv( geneMatrixFile, sep='\t')

if __name__ == "__main__":
    main()

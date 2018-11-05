from __future__ import print_function

########################################################################
# File: jbToDRIMM.py
#  executable: jbToDRIMM.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 10/26/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
from tqdm import *

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
        self.parser = argparse.ArgumentParser(description = ' jbToDRIMM - convert juncBASE step6 AS table to DRIMMSeq table.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -i step6_AS_lenNorm.txt --as_type (e.g. cassette)')
        # Add args
        self.parser.add_argument('-i', "--input_jbTable", action = 'store', required=True, help='juncbase table!!')
        self.parser.add_argument("--as_type", action = 'store', default = False, required=False, help='Optional: If you only want to get an event type.')
        self.parser.add_argument("--keep_jcn_only", action = 'store_false', required=False, default = True, help='Filter JCN_ONLY JB entries.')
        self.parser.add_argument("--quiet", action = 'store_false', required=False, default = True, help='Do not display progress')
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

#######
# MAINE
#######

def main():
    '''
    not the state.
    '''
    myCommandLine = CommandLine()
    jbTable       = myCommandLine.args['input_jbTable']
    asType        = myCommandLine.args['as_type']
    jcnOnly       = myCommandLine.args['keep_jcn_only']

    if asType == False:
        filterFlag = False
    global verbose
    verbose = myCommandLine.args['quiet']

    #convert juncbase to drimm
    '''
    format follow:
    feature_id gene_id samples....
    '''

    count = 0
    if verbose:
        try:
            print("Reading JuncBASE table %s"  % jbTable, file=sys.stderr)
            with open(jbTable) as lines:
                for line in lines:               
                    count += 1
        except:
            print("Cannot read JuncBASE table %s" % jbTable, file=sys.stderr)
    else:
        pass


    jbAmended = open("ammended_jb_table.tsv",'w')
    fileDict = dict()
    with open(jbTable) as lines:

        header = next(lines)
        header = header.rstrip().split("\t")
        eventUUID = 0 
        for line in tqdm(lines, total=count, desc="Converting juncbase entries to DRIMMSeq format.") if verbose else lines:
    
            cols       =  line.rstrip().split("\t")
            gene       = cols[2]
            asType     = cols[1]
            novelty    = cols[0]
            uuidCoords = cols[5]

            if len(gene)<1:
                gene = uuidCoords

            if (filterFlag and cols[1] != asType) or (jcnOnly and "jcn_only" in asType):
                continue
               
            vals = cols[11:]
            inclusions = [x.split(";")[0] for x in vals]
            exclusions = [x.split(";")[1] for x in vals]
            eventUUID += 1

            if asType not in fileDict:
                fileDict[asType] = open("jb_converted_drim_table_%s.tsv" % asType,'w')
                print("feature_id","gene_id","\t".join(header[11:]), sep="\t",file=fileDict[asType])

            print("%s_inclusion_%s" % (uuidCoords, eventUUID), "%s;%s;%s;%s;%s" % (gene,uuidCoords,asType,novelty,eventUUID), "\t".join(inclusions),sep="\t",file=fileDict[asType])
            print("%s_exclusion_%s" % (uuidCoords, eventUUID), "%s;%s;%s;%s;%s" % (gene,uuidCoords,asType,novelty,eventUUID), "\t".join(exclusions),sep="\t", file=fileDict[asType])
    
            print("%s;%s;%s;%s;%s" % (gene,uuidCoords,asType,novelty,eventUUID), line.rstrip(), sep="\t", file=jbAmended)
            

    for k,v in fileDict.items():
        v.close()
    jbAmended.close() 
if __name__ == "__main__":
    main()

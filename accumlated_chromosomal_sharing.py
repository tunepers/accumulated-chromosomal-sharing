#!/usr/bin/python

import pdb
from datetime import date
import pandas as pd
from bx.intervals.cluster import ClusterTree
from bx.intervals.intersection import Interval
from bx.intervals.intersection import IntervalTree
import argparse


parser = argparse.ArgumentParser(description='Get coverage')
parser.add_argument('input_overlap_file', metavar='Input PLINK overlap file', type=str, help='Input file from PLINK for a given chromosome')
parser.add_argument('input_fam_file', metavar='Input PLINK .fam file', type=str, help='Input file PLINK .fam file ')
args = parser.parse_args()


####### Read overlap file
collection_dict = {}
with open(args.input_overlap_file,'r') as infile:
    for i,line in enumerate(infile.readlines()[1:]):
        words = line.strip().split()
        if not 'BP1' in line:

            if not words[1] in collection_dict:
                collection_dict[words[1]] = {} 
            collection_dict[words[1]][words[3]] = (int(words[6]),int(words[7]),i) # "%s-%s-%s"%(int(words[6]),int(words[7]),i)
        
            if not words[3] in collection_dict:
                collection_dict[words[3]] = {} 
            collection_dict[words[3]][words[1]] =  (int(words[6]),int(words[7]),i) # "%s-%s-%s"%(int(words[6]),int(words[7]),i) 



####### Function to sum of region (in terms of bps)

def get_region_bps_sum(tree):
    bps_sum = 0
    for start, end, name in tree.getregions():
        bps_sum += end-start
    return bps_sum


####### Identify overlaps
ind = pd.read_csv(args.input_fam_file,sep="\t")
results_df = pd.DataFrame(index=list(ind.IID), columns=range(1,len(ind)+1))
for i,IID1 in enumerate(ind.IID):
    
    # Status
    if ( i % 100 == 0):
        print "%s out of %s"%(i+100, len(ind.IID))

    tree = ClusterTree(0,0)
    running_sum = 0
    
    # All entries for IID and IID2
    for j,IID2 in enumerate(ind.IID):

        if IID1 in collection_dict and IID2 in collection_dict[IID1]:
            sta, end, i = collection_dict[IID1][IID2]
            tree.insert(sta, end, i)
            running_sum = get_region_bps_sum(tree)

        results_df.at[IID1,j+1] = running_sum

####### Write to file
results_df.to_csv("%s_overlap.txt"%infilename, index=True, quoting=0, doublequote = False, sep='\t')

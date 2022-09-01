#!/usr/bin/python3

'''
This script will load BLAST results in the tabular format and identify reciprocal best hit genes.

'''

import sys
from operator import itemgetter, attrgetter
from pprint import pprint

# load two species into two lists
spec1 = []
spec2 = []
spec1n = sys.argv[2]
spec2n = sys.argv[3]
outfn = sys.argv[4]
for line in open(sys.argv[1],'r'):
    if line.find(spec1n) == 0:
        spec1.append(line.strip())
    if line.find(spec2n) == 0:
        spec2.append(line.strip())

def getgenename(input):
    # this is a parser that convert species specific isoform names
    # to gene names.
    #
    if input.find('ARATH')==0:
        return 'ARATH|'+input.split('|')[1].split('.')[0]
    elif input.find('GLYMA')==0:
        return 'GLYMA|'+'.'.join(input.split('|')[1].split('.')[0:2])
    else:
        return 'UnKnownSpecies'

def getbesthit(speclist):
    allhits = {} # this holds all BLAST hits for each gene
    for eg in speclist:
        tmp = eg.split('\t')
        queryname = getgenename(tmp[0])
        targetname = getgenename(tmp[1])
        eval = float(tmp[10])
        if not queryname in allhits:
            allhits[queryname]=[]
        allhits[queryname].append([targetname,eval])
    # sort allhits and get best ones.
    besthits= {}
    for eg in allhits:
        tmp = allhits[eg]
        besthita=sorted(tmp, key=itemgetter(1))
        besthitg=besthita[0]
        besthits[eg] = besthitg
        #pprint(besthita)
        #pprint(besthitg)
        #print eg,besthitg

    return besthits

spec1bh = getbesthit(spec1)
spec2bh = getbesthit(spec2)

#pprint(spec1bh)
#pprint(spec2bh)

outf = open(outfn,'w')
for eg in spec1bh:
    bhg = spec1bh[eg][0] # name of best hit
    bhe = spec1bh[eg][1] # e value for best hit
    if not bhg in spec2bh:
        continue # skip this if best hit is not even in other species list
    else:
        rbhg = spec2bh[bhg][0] # find best hit gene in species 2
        rbhe = spec2bh[bhg][1] # and associated e value.
        if rbhg == eg:
            outf.write(bhg+','+rbhg+','+`bhe`+','+`rbhe`+'\n')







#!/usr/bin/python3
import sys

TOPN = int(sys.argv[1])
print "=== Top "+str(TOPN)+" One Way Best Hit results ==="

import numpy as np
inFile = sys.argv[2]
outFile = inFile[0:-4]+"_Top-"+str(TOPN)+"_OneWayBestHit.txt"
print "- Output file    : ", 
print outFile

with open(inFile, "r") as f:
	blastList=[]

	for line in f:
		bList=[]
		line=line.strip("\n\r").split('\t')
		qr= '.'.join(line[0].split(".")[0:-1])
		db= '.'.join(line[1].split(".")[0:-1])
		evl=float(line[10])
		bit=float(line[11])
		lth=int(line[3])
		idn=float(line[2])
		bList=[qr, db, evl, bit, lth, idn]
		blastList.append(bList)

print "- blastList All  : ", 		
print len(blastList)

import itertools
blastList.sort()
blastListUnq=list(blastList for blastList,_ in itertools.groupby(blastList))

print "- blastList Uniq : ", 
print len(blastListUnq)

blastListUnqArray=np.array(blastListUnq, dtype=object)

print "- length of Query: ", 
queries=list( set( np.array(blastListUnq)[:,0] ) )
print len(queries)

headerList=[["Query", "DB", "E-val", "Bits", "length", "PerctIdentity"]]
blastArrayTop5s = np.array(headerList)

#print blastArrayTop5s
for q in queries:
	blastQarray = blastListUnqArray[np.where(blastListUnqArray[:,0] == q)].tolist()
	blastArraySort= np.array(sorted(sorted(sorted(sorted(sorted(blastQarray,key=lambda e:e[5],reverse=True),key=lambda e:e[4],reverse=True),key=lambda e:e[3],reverse=True),key=lambda e:e[2],reverse=False),key=lambda e:e[0],reverse=False))[0:TOPN]
	blastArrayTop5s=np.concatenate((blastArrayTop5s, blastArraySort),axis=0)

print "- Num of results : ", 
print len(blastArrayTop5s)-1   #exclude the header line for counting lines.

import pandas as pd 
df = pd.DataFrame(blastArrayTop5s)
df.to_csv(outFile, header=None, index=False, sep="\t")


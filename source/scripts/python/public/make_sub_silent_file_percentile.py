#!/usr/bin/python

import string
from sys import argv,stderr
from os import popen
import sys
import os

def Help():
    print '\n'+'-'*75
    print 'usage: %s <silent-file> <out-file> <score-index> <number-to-choose>\n'%argv[0]
    print 'negative score-index means choose lowest'
    print 'score-index is as in sort; 1=total_score, etc'
    print 'score-index = 0 will choose a random subset'
    print '-'*75
    print '\n'
    sys.exit(0)
    return

if len(argv) != 5:
    Help()

in_file = argv[1]

if not os.path.isfile( in_file ):
	print "Cannot find inputfile for make_sub_silent_file_percentile.py: ", in_file 
	sys.exit(1)


out_file = argv[2]
index = int(argv[3])
choose = float(argv[4])

#########################

lines = map(string.split,
            popen('grep "SCORE" '+in_file+ " | grep -v \" score \" ").readlines()[1:]) ## skip header

N = len(lines)

choose = choose * N
choose = int(choose)
print N, choose

if N <= choose:
  choose = N-1

scores = []
for line in lines:
    if len(line) >= 15:
        try:
          the_score = float(line[abs(index)])
        except:
          the_score = 100000000.0
        scores.append( the_score )
scores.sort()
#  stderr.write('score_index: %d min: %f 1st: %f 10th: %f 90th: %f 99th: %f max: %f\n'\
#                   %(index,
#                     scores[0],scores[N/100],scores[N/10],
#                     scores[(90*N)/100],scores[(99*N)/100],scores[N-1]))

if index<0:
    multiplier = -1
    index = -1*index
else:
    multiplier = 1
    scores.reverse()
threshold = multiplier * scores[choose]
stderr.write('threshold: %f\n'%threshold)
    

    
data = open(in_file,'r')
out = open(out_file,'w')

out.write(data.readline())
out.write(data.readline())
line = data.readline()
counter = 0
nwrites = 0
while line:
    if line[:5] != 'SCORE':
        line = data.readline()
        continue

    counter = counter+1
    if not counter%1000:stderr.write('%d of %d\n'%(counter,N))

    if nwrites >= choose: break
    
    l = string.split(line)
    if len(line)<15 or l[-1] == 'description':
        line = data.readline()
        continue
    write = 0
   
    try:
      the_score = float(l[index])  
    except:
      the_score = 100000000.0
     
    if multiplier > 0:
      if multiplier * the_score >= threshold:
        write = 1
    else:
      if multiplier * the_score >= threshold:
        write = 1
 

    if write:
        out.write(line)
        nwrites = nwrites + 1
    line = data.readline()
    while line and line[:5] != 'SCORE':
        if write:
            out.write(line)
        line = data.readline()
data.close()
out.close()

#!/usr/bin/env python2.7

import sys
import re

def main(argv):
	if(len(argv) == 0 or argv[0] == "-h"):
		print "usage: combine_ContactMaps.py file1 file2 ... > combined_file"
		return
	matrix=[]
	models=0
	first_file = open(argv[0], 'r')
	for line in first_file:
		line=line.rstrip()
		if re.match("#", line) :
			words = line.split()
			models += int (words[-1]) 
		elif line:
			matrix.append(line.split('\t'))
	
	for arg in argv[1:]:
		file = open(arg, 'r')
		nrow=0
		for line in file:
			line=line.rstrip()
			if re.match("#", line) :
                   		words = line.split()
                        	models += int (words[-1])
				continue
			elif not line :
				continue
			if nrow== 0:
				nrow+=1
				continue
			ncol=0
			for col in line.split('\t'):
				if ncol==0:
					ncol += 1
					continue
			#	print "col" +ncol + ", row" +nrow+": "+col
			#	print matrix[nrow][ncol]
				matrix[nrow][ncol] = int(matrix[nrow][ncol]) + int(col)
				ncol += 1
			nrow +=1
	print '# Number of Models: {}'.format(models)
	print 
	for row in matrix:
		output_string= str(row[0])
		for col in row[1:]:
			 output_string += '\t' + str(col)  	
		print output_string

if __name__ == "__main__":
    main(sys.argv[1:])

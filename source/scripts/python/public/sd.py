#!/usr/bin/env python

#######################################################################################################################################################
def syntax():
	print "sd 0.23 (SilentData) Extracts data columns from silent files. (c) Mike Tyka 2005-2011 mike.tyka@gmail.com "
	print "Syntax: sd file1 [file2 file3 .. ] -- <columnname1> [  <columnname2> ... ] "
	print "Examples: "
	print " sd file.out                      ## extract rms and score column (by default)"
	print " sd file.out gdtmm rms score      ## extract gdtmm, rms and score"
	print " sd file.out score                ## extract score only"
	print "\nAdvanced:"
	print " sd file1.out file2.out -         ## extract rms and score from file1.out file2.out (sequentially) "
	print "                                  ## Note the '-' character MUST be used if multiple files are specified"
	print " sd *.out -                       ## extract rms and score from all silent files in directory"
	print " sd *.out *.out.gz *.out.bz2 -    ## extract rms and score from all files including gzipped and bziped ones (sd has built in on-the-fly decompression)"
	print " sd *.out - gdtmm rms description ## extract gdtmm, rms and description"
	print " sd *.out - score __file__        ## extract score and the filename"
	print "\nNinja:"
	print " sd *.out - 1.0 score -0.2 fa_atr 0.3 fa_sol      ## extract from all outfiles the sum of (1*score - 0.2*fa_atr + 0.3* fa_sol )"
	print " sd *.out - rms 1.0 score -0.2 fa_atr 0.3 fa_sol  ## extract from all outfiles rms and the sum of (1*score - 0.2*fa_atr + 0.3* fa_sol )"
#######################################################################################################################################################

import sys
import os
import string
import re
import gzip
import bz2

non_numerical_cols = [
"aln_id",
"user_tag",
"description",
"usid",
"husid",
"__file__"
];


if len(sys.argv) < 2:
	syntax()
	sys.exit(1)

## Analyse options

showfilename = False
multifilemode = False
summode = False
normalize_sums = False
normalize_all = False
offset_sums = False
offset_all = False

for arg in sys.argv[1:]:
	if arg == "-":
		multifilemode = True
		continue


columnbreak = False

columns = []
files = []
sums = []


next_sum_type = ""

for arg in sys.argv[1:]:
	if arg == "-o":
		offset_all= True
		continue
	if arg == "-os":
		offset_sums = True
		continue
	if arg == "-n":
		normalize_all = True
		continue
	if arg == "-ns":
		normalize_sums = True
		continue
	if arg == "-f":
		showfilename = True
		continue
	if arg == "-":
		columnbreak = True
		continue

	if not columnbreak:
		files.append( arg )
		if not multifilemode:
			columnbreak = True
		continue

	if arg[0] == '-' or arg[0] == '+' or arg[0] == '.' or (arg[0] >= '0' and arg[0] <= '9'):
		next_sum_type = arg
		summode = True
		continue

	column_name = arg
	sum_type = next_sum_type
	columns.append( ( column_name, sum_type ) )

	next_sum_type = ""

if len(columns) == 0: columns = [('rms',''), ('score','') ]

final_data = []


## now iterate through files and print columns etc..
for filename in files:

	## gzip/bzip support
	lines = []
	if( filename[-3:] == ".gz" ): ## gzipped
		file = gzip.open(filename,"r")
		lines = string.split(file.read(),"\n")
	elif( filename[-4:] == ".bz2" ): ## gzipped
		lines = bz2.BZ2File( filename,"r").readlines()
	else: ## plain text
		file = open(filename,"r")
		lines = string.split(file.read(),"\n")


	donetitles = False
	names = []
	for arg,sumtype in columns:
		names.append((arg,-1))

	for (i,(n,j)) in enumerate(names):
		names[i] = (string.lower(n), j)

	if len(names)<=0:
		sys.stderr.write("Warning: no columns names supplied: defaulting to rms vs score\n")
		names.append( ("rms", -1) )
		names.append( ("score", -1) )


	for number,l in enumerate(lines):
		token = string.split(l)
		if len(token) <= 0: continue
		## SCORE LINE HEADERS
		if not donetitles:
		 if (token[0] == "SCORE ") or (token[0] == "SCORE:") :
			## look for names

			for ti,t in enumerate(token):
				for (i,(n,j)) in enumerate(names):
					if string.lower(t) == n:	## if names match
						(existing_n, existing_ti) = names[i];
						if( existing_ti < 0):
							names[i] = (n,ti)	## remember index



			Error = False
			for i,(n,j) in enumerate(names):
				if n == "__file__": ## special column
					names[i] = (n, -1 )
				elif j < 0:
					sys.stderr.write("Error: Cannot find column named '%s' in '%s' \n"%(n,filename) )
					Error = True

			if Error: sys.exit(1)
			donetitles = True
			continue


		## SCORE LINE
		if token[0] == "SCORE:":
			token = string.split(l)
			brokenline = False
			for (n,i) in names:
				if i >= len( token ):
					brokenline = True
					continue
			if brokenline : continue

			## Check columns for correctness. Allow only special columns to have
			# non-numerical values ( see array at the top of this code file)
			for (n,i) in names:
				## Figure out if the tag is an official non-numeric one
				is_nonnumtag = False
				for nonnumtags in non_numerical_cols:
					if n == nonnumtags:
						is_nonnumtag = True
						break
				if is_nonnumtag: break

				##print n,i,token[i]
				## if we're here, the tag is considered a numeric one and thus must pass the float conversion test
				try:
					t = float( token[i] )
				except:
				  ## if it does not then ignore this line - we assume it's corrupted
					brokenline = True;

			## Ignore broken lines!
			if brokenline : continue

			## Ok, if we're here, the line is clean and we can go ahead to extract the data
			newdataline = dict();

			# now actuall extract the data into an array
			for (n,i) in names:
				if n == "__file__":
					## record the filename as "data"
					newdataline[ n ] = filename
				else:
					## extract the data itself
					newdataline[ n ] = token[i]
				#sys.stdout.write("%10s "%token[i])
			#sys.stdout.write("\n")


			## and append that array to the overall data stack
			final_data.append( newdataline )


################################################

# Now print combined data

combined_data = []

for l in final_data:
	sum = 0
	sumactive = False

	new_combined_data_line = []

	for (n,i) in columns:

		if i == '':
			if sumactive:
				#sys.stdout.write("%10s "%sum )
				new_combined_data_line.append( ("", sum  ) )
				sum = 0
				sumactive = False

			#sys.stdout.write("%10s "%l[ n ])
			new_combined_data_line.append( ( n, l[ n ] ) )
		else:
			sum += float(l[ n ]) * float(i)
			sumactive = True

	if sumactive:
			#sys.stdout.write("%10s "%sum )
			new_combined_data_line.append( ("", sum ) )
			sum = 0
			sumactive = False

	combined_data.append( new_combined_data_line )


range_data = dict()

for l in combined_data:
	for (n,i) in l:
		if range_data.get( n ): range_data[n].append( i )
		else: range_data[n] = [ i ];


perc_high = dict()
perc_low  = dict()


for key in range_data.keys():
	range_data[key].sort()
	perc_high[key] = range_data[key][ int( len( range_data[key] ) * 0.95)  ]
	perc_low[key]  = range_data[key][ int( len( range_data[key] ) * 0.05)  ]


for l in combined_data:
	for (n,i) in l:
		try:
			value =  float(i)
			if normalize_all or ( normalize_sums and n == '' ):
				value = ( float(i) - perc_low[ n ] ) / ( perc_high[ n ] - perc_low[ n ] )
			elif offset_all or ( offset_sums and n == '' ):
				value = ( float(i) - perc_low[ n ] )
			sys.stdout.write("%10.3f "%value )
		except:
			sys.stdout.write("%10s "%i )

	sys.stdout.write("\n")



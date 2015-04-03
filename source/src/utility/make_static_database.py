import sys;
import os
import os.path
import re

output = "static_database.hh";
path_prefix = "/Users/mtyka/minirosetta_database/"

outfile = open( output, "w" )

outfile.write('#ifndef STATIC_DATABASE_HH\n')
outfile.write('#define STATIC_DATABASE_HH\n')

outfile.write('const char* static_database[ ][2] = {\n')

# some regulat expressions we'll need

escaped_chars = re.compile('(["\\%])');


for file in sys.argv[1:] :

	if( file[0:2] == "./" ): file = file[2:]

	file = re.sub( '//', '/', file )


	
	print "Reading: ", file


	try:
		input = open( path_prefix + "/" +  file, "r" );
	except: 
		outfile.write('{"'+file+'",""},\n')
		continue

	outfile.write('{\n')
	outfile.write('"'+file+'",\n"\\\n')
	inputdata = input.read().split("\n")
	for line in inputdata:
		if line == "": continue
		#prefix these characters with backslashes
		line = re.sub( '\\\\', '\\\\\\\\', line )  
		line = re.sub( '"', '\\"', line )  
		line = re.sub( '%', '\\%', line )  
		outfile.write( line + '\\n\\\n' )
	outfile.write( '"},\n' )

outfile.write('{"",""}};\n')
outfile.write('const int static_database_size = %d;\n\n'%(len(sys.argv[1:])) )
outfile.write('#endif\n\n')
outfile.close()

#!/usr/bin/env python

from optparse import OptionParser
import os.path
import sets
import sys

"""
This script is meant to facilitate updating the boost c++ header-only libraries
in mini.  Given a boost library and its primary header filenames/directory,
it will run through and find the other headers necessary for operation.  Run
without command line arguments to obtain valid script options.  At the end of
a run, the script will indicate how to use rsync to copy the headers.  Desired
boost libraries and their primary header filenames/directories are stored
below in the global variable BOOST.  If you wish to add a boost library please
modify the BOOST variable.
"""

mVersion = 0.001

###########
# methods #
###########
def boost_path( filename ):
	return os.path.join( "boost", filename )

def boost_path( top_dir, filename ):
	return os.path.join( top_dir, "boost", filename )

def boost_includes( filename ):
	incl = []
	f = open( filename )
	for line in f:
		line = line.strip()
		if len( line ) > 0 and line[ 0 ] == '#':
			# spacing can be somewhat variable, incrementally
			# trim down the line to get the desired result
			line = line[ 1: ].strip()
			if line[ 0:7 ] == "include":
				line = line[ 7: ].strip().split()[ 0 ] # take filename only
				line = line.lstrip( '<' )
				line = line.rstrip( '>' )
				line = line.lstrip( '"' )
				line = line.rstrip( '"' )
				if line[ 0:6 ] == "boost/":
					incl.append( line )
	f.close()
	return incl

###########################################
# Components for desired boost libraries. #
###########################################
class BoostLib:
	def __init__( self, name, components=[], recursive=True ):
		"""Initializes with name of boost library and components (can be primary header
		   files or directories).

		   If recursive is True, then all files will be searched for includes.
		   If recursive is false, then only files/directories indicated will be
		   output.
		"""
		self.name = name
		self.components = components
		self.recursive = recursive 

	def check( self, top_dir ):
		"""Check to see if primary header files/directories exist."""
		for i in self.components:
			filename = boost_path( top_dir, i )
			if not os.path.exists( filename ):	
				sys.stderr.write( "ERROR: missing %s\n" % ( filename ) )
				return False
		return True

	def find_components( self, top_dir ):
		"""Find all headers required for operation."""
		unprocessed = sets.Set( [ boost_path( top_dir, i ) for i in self.components ] )
		processed = sets.Set()

		while len( unprocessed ) > 0:
			component = unprocessed.pop()
			if os.path.isfile( component ):
				processed.add( component )
				if self.recursive: # find all includes within file
					new_incl = boost_includes( component )
					for i in new_incl:
						ii = os.path.join( top_dir, i )
						if ii not in processed:
							unprocessed.add( ii )
			elif os.path.isdir( component ):
				# add all files inside directory to unprocessed
				for root, dirs, files in os.walk( component ):
					for f in files:
						ii = os.path.join( root, f )
						if ii not in processed:
							unprocessed.add( ii )
			elif os.path.islink( component ): 
				# in principal links can be handled fine
				# as is, but we shouldn't encounter any
				# so fail-fast
				sys.stderr.write( "ERROR: encountered symbolic link, handling not formalized: %s\n" % ( component ) )
				sys.exit( 1 )
		return processed # the set of all includes

## To add a boost library, add a BoostLib object below to BOOST with the
## correct primary header file(s) and directories, if applicable.  If only
## the indicated header file(s) and directories should be copied then
## then set recursive=False when initializing the BoostLib object.
##
## Current as of boost 1.55 integration (feb 2014)
BOOST = [
        BoostLib( "accumulators", [ "accumulators/", "accumulators/statistics/" ] ),
        BoostLib( "algorithm", [ "algorithm/string.hpp", "algorithm/string" ] ),
        BoostLib( "archive", ["archive/"] ),
        BoostLib( "array", ["array.hpp"] ),
        BoostLib( "assign", [ "assign/", "assign/std/", "assign.hpp" ] ),
        BoostLib( "bind", [ "bind.hpp" ] ),
        BoostLib( "config", [ "config.hpp", "config/", "config/compiler/","config/stdlib/","config/platform/"  ] ),
        BoostLib( "cstdint", [ "cstdint.hpp" ] ),
        BoostLib( "date_time", [ "date_time/posix_time/posix_time_types.hpp" ] ),
        BoostLib( "detail", [ "detail/atomic_count.hpp" ] ),
        BoostLib( "dynamic_bitset", [ "dynamic_bitset/dynamic_bitset.hpp" ] ),
        BoostLib( "foreach", [ "foreach.hpp" ] ),
        BoostLib( "format", [ "format.hpp" ] ),
        BoostLib( "function", [ "function.hpp" ] ),
        BoostLib( "functional", [ "functional/factory.hpp", "functional/hash.hpp" ] ),
        BoostLib( "graph", [ "graph/" ] ),
        BoostLib( "io", [ "io/" ] ),
        BoostLib( "lexical_cast", [ "lexical_cast.hpp" ] ),
        BoostLib( "math", [ "math/constants", "math/distributions", "math/distributions.hpp", "math/special_functions" ] ),
        BoostLib( "mpi", [ "mpi.hpp", "mpi/" ] ),
        BoostLib( "mpl", [ "mpl/" ] ),
        BoostLib( "noncopyable", [ "noncopyable.hpp" ] ),
        BoostLib( "numeric", [ "numeric/ublas" ] ),
        BoostLib( "optional", [ "optional.hpp" ] ),
        BoostLib( "pool", [ "pool/", "pool/detail/" ] ), 
        BoostLib( "preprocessor", [ "preprocessor/preprocessor.hpp" ] ),
        BoostLib( "progress", [ "progress.hpp" ] ),
        BoostLib( "python", [ "python.hpp", "python/" ] ),
        BoostLib( "scoped_ptr", [ "scoped_ptr.hpp" ] ),
        BoostLib( "serialization", [ "serialization/", "serialization/detail/" ] ),
        BoostLib( "smart_ptr", [ "smart_ptr/", "smart_ptr.hpp" ] ),
        BoostLib( "spirit", [ "spirit/", "spirit/include/", "spirit/home/phoenix/bind/bind_function.hpp" ] ),
        BoostLib( "thread", [ "thread.hpp", "thread" ] ),
        BoostLib( "timer", [ "timer.hpp" ] ),
        BoostLib( "tokenizer", [ "tokenizer.hpp" ] ),
        BoostLib( "tuple", [ "tuple" ] ),
        BoostLib( "type_traits", [ "type_traits/", "type_traits.hpp" ] ),
        BoostLib( "unordered", [ "unordered/", "unordered_map.hpp", "unordered_set.hpp" ] ),
        BoostLib( "utility", [ "utility/", "utility.hpp" ] ),
        BoostLib( "uuid", [ "uuid" ] ),
        BoostLib( "variant", [ "variant.hpp" ] ),
        BoostLib( "version", [ "version.hpp" ] ),
        BoostLib( "weak_ptr", [ "weak_ptr.hpp" ] ),
        BoostLib( "xpressive", [ "xpressive" ] ),
]


########
# main #
########
def main( options, args ):
	# process options
	if not options.boost_out:
		sys.stdout.write( "* No output directory specified, entering checking mode to verify library list.\n" )

	# run through all desired libraries and check first
	for lib in BOOST:
		sys.stdout.write( "* checking: %s\n" % lib.name )
		if not lib.check( options.boost_in ):
			sys.exit( 1 )

	# now do actual copying
	if options.boost_out:
		# gather list of all necessary headers
		required = sets.Set()
		for lib in BOOST:
			sys.stdout.write( "* processing: %s\n" % lib.name )
			required.update( lib.find_components( options.boost_in ) )
		# write list of necessary headers 
		out = open( options.boost_out, 'w' )
		for i in required:
			relpath = i.replace( options.boost_in + os.sep, '', 1 )
			out.write( relpath + '\n' )
		out.close()

		# notify user how to copy headers
		sys.stdout.write( "* Required boost includes written to: %s\n" % ( options.boost_out ) )
		sys.stdout.write( "* Now use rsync to copy the files (must use '-ar' flags),\n" )
		sys.stdout.write( "*    rsync -ar --files-from=%s %s NEW_DIR\n" % ( options.boost_out, options.boost_in ) )

	sys.stdout.write( "* Done.\n" )

##############
# invocation # 
##############

if __name__ == '__main__':
	# setup options parser
	usage = "usage : %prog [OPTIONS]"
	parser = OptionParser( usage, version=mVersion )
	parser.add_option( "-i", "--boost_in", action="store", type="string", dest="boost_in", default=None, help="take boost headers from this *top-level* directory, e.g. /usr/include or boost_1_38_0 (REQUIRED)" )
	parser.add_option( "-o", "--out", action= "store", type="string", dest="boost_out", default=None, help="write required boost includes to this file; without this option program enters checking mode" )

	# parse arguments
	options, args = parser.parse_args()

	# check for required options
	flag = False
	if not options.boost_in:
		flag = True

	# launch
	if flag:
		parser.print_help()
		sys.exit( 1 )
	else:
		main( options, args )

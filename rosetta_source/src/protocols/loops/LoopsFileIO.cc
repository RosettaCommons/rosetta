// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopsFileIO.cc
/// @brief
/// @author Brian D. Weitzner

// Unit header
#include <protocols/loops/LoopsFileIO.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh> // needed to parse old style loop files
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>




namespace protocols {
namespace loops {

static basic::Tracer tr("protocols.loops.LoopsFileIO");

LoopsFileIO::LoopsFileIO() : utility::pointer::ReferenceCount()
{

}

LoopsFileIO::LoopsFileIO( const LoopsFileIO & ) : utility::pointer::ReferenceCount()
{
}

// destructor
LoopsFileIO::~LoopsFileIO(){}

//////////////////////////////////////////////////////////////////////
std::ostream & operator<< ( std::ostream & os, const LoopsFileIO & /*loops*/ ) {
	/*
	os << "LOOP  begin  end  cut  skip_rate  extended" << std::endl;
	for ( Loops::const_iterator it = loops.begin(), it_end = loops.end();
		 it != it_end; ++it ) {
		os << *it << std::endl;
	}
	*/
	return os;
}


void LoopsFileIO::read_stream_to_END(
	std::istream & is,
	std::string filename /*for error msg */,
	bool strict_looprelax_checks,
	std::string token )
	{
	std::string line;
	int linecount=0;
	int errcount=50; //if we reach 0 we bail!

	loops_.clear();

	while( getline( is, line) ) {
		linecount++;
		utility::vector1< std::string > tokens ( utility::split( line ) );

		SerializedLoop current_loop;
		if( tokens.size() > 0 ) {
			if ( tokens[1].substr(0,3) == "END" ) break;
			if ( tokens[1] == token ) {
				if ( tokens.size() < 3 ) {
					utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Minimum of 3 tokens necessary (begin, end, cutpoint)"  );
				}
				if ( tokens.size() > 6 ) {
					utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Maximum of 6 tokens allowed (LOOP begin end cutpoint skiprate extended)"  );
				}
				current_loop.start = (core::Size) atoi(tokens[2].c_str());
				current_loop.stop = (core::Size) atoi(tokens[3].c_str());
				current_loop.cut = 0;        // default - let LoopRebuild choose cutpoint
				current_loop.skip_rate = 0.0;  // default - never skip
				std::string extend_loop_str;
				bool extend_loop = false;

				if (tokens.size() > 3)
					current_loop.cut = (core::Size) atoi(tokens[4].c_str());
				if (tokens.size() > 4)
					current_loop.skip_rate = atof(tokens[5].c_str());
				if (tokens.size() > 5){
					if( tokens[6] == "X" ){
						tr.Error << "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "[WARNING] DEPRECATED old style extended marker X is used" << std::endl;
						extend_loop = true;
						if ( errcount > 0 ) errcount--;
						else {
							utility_exit_with_message( "too many errors in loop-file " + filename );
						}
					}else{
						int extended_token = atoi(tokens[6].c_str());
						if( extended_token == 0 ) extend_loop = false;
						else                      extend_loop = true;
					}
				}

				current_loop.extended = extend_loop;
				if ( current_loop.start > current_loop.stop || ( current_loop.start==current_loop.stop && strict_looprelax_checks ) ) {
					utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Invalid loop definition (start residue " + ( strict_looprelax_checks ? ">=" : ">" )  + " end residue) - ERROR"  );
				} else {
					loops_.push_back( current_loop );
				}
			} else if ( tokens[1][0] != '#' ) {
				if (tokens.size() >= 2) {
					tr.Error << "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "DEPRECATED r++ style loopfile" << std::endl;

					if ( errcount>0 ) errcount--;
					else {
						utility_exit_with_message( "too many errors in loop-file " + filename );
					}

					current_loop.start = (core::Size) atoi(tokens[1].c_str());
					current_loop.stop   = (core::Size) atoi(tokens[2].c_str());
					current_loop.cut = 0;        // default - let LoopRebuild choose cutpoint
					current_loop.skip_rate = 0.0;  // default - never skip
					bool extend_loop = false;    // default - not extended
					if (tokens.size() > 2)
						current_loop.cut = (core::Size) atoi(tokens[3].c_str());
					if (tokens.size() > 3)
						current_loop.skip_rate = atof(tokens[4].c_str());
					if (tokens.size() > 4){
						if( tokens[5] == "X" ){
							tr.Error << "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "[WARNING] DEPRECATED old style extended marker X is used" << std::endl;
							extend_loop = true;
						} else {
							int extended_token = atoi(tokens[5].c_str());
							if ( extended_token == 0 ) extend_loop = false;
							else                extend_loop = true;
						}
					}
					current_loop.extended = extend_loop;


					if ( current_loop.start > current_loop.stop || ( current_loop.start==current_loop.stop && strict_looprelax_checks ) ) {
						utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Invalid loop definition (start residue " + ( strict_looprelax_checks ? ">=" : ">" ) + "end residue) - ERROR"  );
					}

					loops_.push_back( current_loop );

				} else {
					tr.Warning << "[WARNING] Skipping line '" << line << "'" << std::endl;
				}
			}
 		}
	} //while
}

SerializedLoopList LoopsFileIO::read_loop_file( std::string filename )
{
	std::ifstream infile( filename.c_str() );

	if (!infile.good()) {
		utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + filename + "'" );
	}

	// TODO: The read_stream_to_END stuff should be moved to its own class that will handle the various types of input
	// files
	read_stream_to_END( infile, filename );

	return loops_;
}

SerializedLoopList LoopsFileIO::use_custom_legacy_file_format(
	std::istream & is,
	std::string filename,
	bool strict_looprelax_checks,
	std::string token )
{
	read_stream_to_END( is, filename, strict_looprelax_checks, token );
	return loops_;
}

LoopsFileIO & LoopsFileIO::operator =( LoopsFileIO const & src )
{
		loops_ = src.loops_;
		return *this;
}


} // namespace loops
} // namespace protocols

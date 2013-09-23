// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson


#include <core/types.hh>
#include <devel/init.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

// C++ headers
#include <map>
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


///////////////////////////////////////////////////////////////////////////////

using std::map;
using std::string;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using namespace basic;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int
main( int argc, char* argv [] )
{
	try {

	using core::Real;
	using core::Size;

	// options, random initialization
	devel::init( argc, argv );

	typedef vector1< Sequence > seqlist;
	seqlist seqs( read_fasta_file( option[ in::file::fasta ]()[1] ) );

	Size total( 0 );
	map< char, Size > counts;
	for ( seqlist::iterator it = seqs.begin(), end = seqs.end();
				it != end; ++it
	) {
		for ( Size i = 1; i <= it->length(); ++i ) {
			char const aa_name( (*it)[i] );
			if ( !core::chemical::oneletter_code_specifies_aa( aa_name ) ) continue;
			counts[ aa_name ]++;
			++total;
		}
	}

	utility::io::ozstream out( "stats.txt" );
	Size const width( 12 );
	Size const precision( 6 );

	out 	<< A( width, "aa" )
			<< A( width, "count" )
			<< A( width, "prob" )
			<< std::endl;
	for ( map< char, Size >::const_iterator it = counts.begin(), end = counts.end();
			it != end; ++it
	) {
		Size count( it->second );
		Real prob ( static_cast< Real > (count) / total );
		out 	<< A( width, it->first )
				<< I( width, count )
				<< F( width, precision, prob )
				<< std::endl;
	}
	out.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0
} // int main( int argc, char * argv [] )

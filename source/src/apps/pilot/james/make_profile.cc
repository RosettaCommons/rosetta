// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file make_profile.cc
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

#include <basic/options/keys/out.OptionKeys.gen.hh>
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
	seqlist query_seqs( read_fasta_file( option[ in::file::fasta ]()[1] ) );
	seqlist seq_db    ( read_fasta_file( option[ in::file::fasta ]()[2] ) );

	Size count( 0 );
	//Size max_compare( 1000 );
	NWAligner nw_aligner;
	SWAligner sw_aligner;
	Real gap_open       = -8;
	Real gap_extend     = -1;
	FileName matrix_fn( "BLOSUM62" );
	ScoringSchemeOP ss( new MatrixScoringScheme( gap_open, gap_extend, matrix_fn ) );

	//Real threshold( 0 );
	vector1< Real > scores;
	for ( seqlist::iterator q_it = query_seqs.begin(), q_end = query_seqs.end();
				q_it != q_end; ++q_it
	) {

		SequenceOP query_seq( new Sequence( *q_it ) );
		if ( query_seq->sequence().find( 'X' ) != string::npos ) continue;

		//std::cout << "query_seq = " << (*query_seq) << std::endl;
		for ( seqlist::iterator db_it = seq_db.begin(), db_end = seq_db.end();
					db_it != db_end; ++db_it
		) {
				SequenceOP db_seq( new Sequence( *db_it ) );
				if ( db_seq->sequence().find( 'X' ) != string::npos ) continue;

				//std::cout << "db_seq = " << (*db_seq) << std::endl;
				SequenceAlignment global_align = nw_aligner.align( query_seq, db_seq, ss );
				scores.push_back( global_align.score() );
				//if ( global_align.score() > threshold )
				//	std::cout << "global alignment: " << global_align << std::endl;

				//SequenceAlignment local_align = sw_aligner.align( query_seq, db_seq, ss );
				//std::cout << "local alignment: " << local_align << std::endl;
		} // db sequences

		++count;
	} // query sequences

	Size const width( 12 );
	Size const precision( 6 );
	utility::io::ozstream out( option[ out::file::silent ]() );
	for ( utility::vector1< Real >::iterator it = scores.begin(), end = scores.end();
				it != end; ++it ) {
		out << F( width, precision, *it ) << std::endl;
	}
	out.close();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
using namespace basic;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;

bool
is_locked( std::string const & fn ) {
	if ( utility::file::file_exists( fn ) ) {
		return true;
	} else {
		utility::io::ozstream out( fn );
		std::cout << "locking " << fn << std::endl;
		out << "locked" << std::endl;
		out.close();
		return false;
	}
}


int
main( int argc, char* argv [] )
{
	try {

	using core::Real;
	using core::Size;
	// options, random initialization
	devel::init( argc, argv );

	Real scoring_threshold( 0.9 );

	NWAligner nw_aligner;
	//Real gap_open  ( -10 );
	//Real gap_extend(  -2 );
	FileName matrix_fn( "BLOSUM62" );

	ScoringSchemeOP ss(
		//new MatrixScoringScheme( gap_open, gap_extend, matrix_fn )
		new SimpleScoringScheme()
	);

	typedef vector1< Sequence > seqlist;
	seqlist seqs( read_fasta_file( option[ in::file::fasta ]()[1] ) );

	map< string, Sequence > represent;
	vector1< string > two_letter_codes;
	for ( seqlist::iterator it = seqs.begin(), end = seqs.end();
				it != end; ++it
	) {
		string pdbid = it->id().substr(0,4);
		string two_letter = pdbid.substr(1,2);
		two_letter_codes.push_back( two_letter );
	}

	Size n_done( 0 );
	Size const n_max( 2000 );
	std::cout << "have " << two_letter_codes.size() << " codes total." << std::endl;
	for ( vector1< string >::const_iterator two_letter = two_letter_codes.begin(),
	 			two_letter_end = two_letter_codes.end();
				two_letter != two_letter_end && n_done <= n_max; ++two_letter
	) {

		std::set< Sequence > unique_seqs;
		string mapname = *two_letter + ".map";
		if ( is_locked(mapname) ) continue;

		for ( seqlist::iterator it = seqs.begin(), end = seqs.end();
					it != end; ++it
		) {
			string pdbid = it->id().substr(0,4);
			string chain = it->id().substr(4,1);
			string this_two_letter = pdbid.substr(1,2);

			if ( *two_letter != this_two_letter ) continue;

			PROF_START( SEQUENCE_COMPARISON );
			Sequence longest_seq( *it );

			for ( seqlist::iterator inner = seqs.begin(), inner_end = seqs.end();
						inner != inner_end; ++inner
			) {
				if ( inner->id() == it->id() ) continue; // skip self-alignments

				string inner_pdbid  = inner->id().substr(0,4);
				string inner_fullid = inner->id().substr(0,5);
				if ( pdbid == inner_pdbid ) {
					//std::cout << "aligning " << std::endl << *it << std::endl << *inner << std::endl;
					SequenceOP this_seq ( new Sequence( *it ) );
					SequenceOP inner_seq( new Sequence( *inner ) );
					SequenceAlignment align_ij = nw_aligner.align( this_seq, inner_seq, ss );
					Real percentage_identity
						= static_cast < Real > (align_ij.identities()) / static_cast < Real > (align_ij.length());
					align_ij.score( percentage_identity );

					//std::cout << align_ij.identities() << " / " << align_ij.length() << std::endl;
					//std::cout << align_ij << std::endl << std::endl;

					if ( align_ij.score() >= scoring_threshold ) {
						// stupid hack: sort lexicographically if sequences have the same size!
						if ( inner->length() == longest_seq.length() &&
							 		inner->id() < longest_seq.id()
						) {
							longest_seq = *inner;
						}

						if ( inner->length() > longest_seq.length() ) {
							longest_seq = *inner;
						}
					} // if score >= threshold
				} // if ( pdbid == inner_pdbid )
			} // for inner

			represent[ it->id() ] = longest_seq;
		} // for seqs

		utility::io::ozstream out( mapname );
		out << "pdbid representative" << std::endl;
		for ( map< string, Sequence >::const_iterator it = represent.begin(), end = represent.end();
					it != end; ++it
		) {
			out << it->first << " " << it->second.id() << std::endl;
			unique_seqs.insert( it->second );
		}
		out.close();

		std::string unique_listfile( *two_letter + ".unique" );
		out.open( unique_listfile );
		out << "pdbid representative" << std::endl;
		for ( std::set< Sequence >::const_iterator it = unique_seqs.begin(), end = unique_seqs.end();
					it != end; ++it
		) {
			out << *it << std::endl;
		}
		out.close();


		std::cout << "done with " << mapname << "." << std::endl;
		++n_done;

		PROF_STOP( SEQUENCE_COMPARISON );
		prof_show();
	} // for two_letter_codes

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file score_aln.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <basic/Tracer.fwd.hh>

#include <utility/excn/Exceptions.hh>


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char* argv [] ) {
	try {

	using std::map;
	using std::string;
	using core::Size;
	using core::Real;
	using utility::vector1;
	using utility::file::FileName;
	using namespace basic;
	using namespace core::sequence;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// options, random initialization
	devel::init( argc, argv );

	// configuration information
	//ScoringSchemeFactory ssf;
	//ScoringSchemeOP ss( ssf.get_scoring_scheme( "Simple" ) );
	ScoringSchemeOP ss( new SimpleScoringScheme( 120, 0, -100, -10 ) );

	SequenceOP remap_seq1(
		read_fasta_file( option[ in::file::fasta ]()[1] )[1]
	);
	//SequenceOP remap_seq2(
	//	read_fasta_file( option[ in::file::fasta ]()[1] )[2]
	//);

	std::string output_fn( option[ out::file::alignment ]() );
	utility::io::ozstream output( output_fn );

	vector1< SequenceAlignment > alns;
	vector1< std::string > align_fns = option[ in::file::alignment ]();

	SWAligner local;
	NWAligner global;

	//Real go = -100;
	//Real ge = -10;
	//ss->gap_open( go );
	//ss->gap_extend( ge );
	//SequenceAlignment al = local.align( remap_seq1, remap_seq2, ss );
	//std::cout << "(" << go << "," << ge << ")" << std::endl;
	//std::cout << al << std::endl;
	//std::cout << "--" << std::endl;
	//std::exit(0);

	typedef vector1< string >::const_iterator aln_iter;
	for ( aln_iter aln_it = align_fns.begin(), aln_end = align_fns.end();
				aln_it != aln_end; ++aln_it
	) {
		alns = core::sequence::read_aln( option[ cm::aln_format ](), *aln_it );

		using namespace core::sequence;
		typedef vector1< SequenceAlignment >::iterator iter;
		for ( iter it = alns.begin(), end = alns.end(); it != end; ++it ) {
			std::cout << "original alignment: " << std::endl << *it << std::endl;

			SequenceOP seq1_copy = it->sequence(1)->clone();
			SequenceOP seq2_copy = it->sequence(2)->clone();
			seq1_copy->sequence( seq1_copy->ungapped_sequence() );
			seq2_copy->sequence( seq2_copy->ungapped_sequence() );

			SequenceAlignment aln2 = local.align( it->sequence(2), remap_seq1, ss );
			//std::cout << aln2 << std::endl;

			SequenceMapping map1 = it->sequence_mapping ( 1, 2 );
			SequenceMapping map2 = aln2.sequence_mapping( 1, 2 );
			SequenceMapping map3 = transitive_map( map1, map2 );
			SequenceAlignment test_aln = mapping_to_alignment(
				map1, seq1_copy, seq2_copy
			);
			std::cout << test_aln << std::endl;
			std::cout << "done!" << std::endl;
			std::exit(0);

			SequenceAlignment new_aln = mapping_to_alignment(
				map3, seq1_copy, seq2_copy
			);

			std::cout << new_aln << std::endl;
			output << new_aln;
			output << "--" << std::endl;

		} // for alns
	} // for align_fns
	output.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

} // int main( int argc, char * argv [] )

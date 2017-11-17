// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file cs_align.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>

#include <protocols/comparative_modeling/AlignmentSet.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/CompositeScoringScheme.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/ChemicalShiftSequence.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

///////////////////////////////////////////////////////////////////////////////

using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using protocols::comparative_modeling::AlignmentSet;
using namespace basic;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;

void print_matrix_to_file(
	std::string const & fn,
	core::sequence::SequenceOP seq1,
	core::sequence::SequenceOP seq2,
	core::sequence::ScoringSchemeOP ss
) {
	std::cout << "printing matrix to " << fn << std::endl;
	utility::io::ozstream out(fn);
	using core::Size;
	for ( Size jj = 1; jj <= seq2->length(); ++jj ) {
		out << " res" << jj;
	}
	out << std::endl;

	for ( Size ii = 1; ii <= seq1->length(); ++ii ) {
		out << "res" << ii;
		for ( Size jj = 1; jj <= seq2->length(); ++jj ) {
			out << " " << F( 8, 3, ss->score(seq1,seq2,ii,jj) );
		}
		out << std::endl;
	}
	out.close();
}

int
main( int argc, char* argv [] ) {
	try {

	using core::Size;
	using core::Real;
	using utility::vector1;
	using core::sequence::SequenceAlignment;

	// options, random initialization
	devel::init( argc, argv );

	Real const gap_open  ( option[ cm::min_gap_open   ]() );
	Real const gap_extend( option[ cm::min_gap_extend ]() );
	std::string const scoring_scheme_type( option[ cm::seq_score ]()[1] );
	ScoringSchemeFactory ssf;
	ScoringSchemeOP ss( ssf.get_scoring_scheme( scoring_scheme_type ) );
	ss->gap_open  ( gap_open   );
	ss->gap_extend( gap_extend );

	std::string output_fn = option[ out::file::alignment ]();
	utility::io::ozstream output( output_fn );

	vector1< FileName > const pssm_fns( option[ in::file::pssm ]() );
	ChemicalShiftSequenceOP seq1( new ChemicalShiftSequence( pssm_fns[1] ) );
	seq1->id( FileName( seq1->id() ).base() );

	std::cout << "starting to generate alignments" << std::endl;
	std::cout << seq1->to_string() << std::endl;

	AlignmentSet align_set;
	vector1< SequenceAlignment > aligns;
	CompositeScoringSchemeOP css(
		new CompositeScoringScheme()
	);
	css->add_scoring_scheme( ss );
	if ( option[ james::debug ]() ) {
		css->add_scoring_scheme( new MatrixScoringScheme(
			gap_open, gap_extend,
			FileName("/work/tex/minirosetta_database/sequence/BLOSUM62") )
		);
	}
	css->gap_open  ( gap_open   );
	css->gap_extend( gap_extend );

	for ( Size jj = 2; jj <= pssm_fns.size(); ++jj ) {
		ChemicalShiftSequenceOP seq2( new ChemicalShiftSequence(pssm_fns[jj]) );
		seq2->id( FileName( seq2->id() ).base() );
		std::cout << seq2->to_string() << std::endl;
		string out_fn( seq2->id() + ".dat" );
		print_matrix_to_file(out_fn,seq1,seq2,css);

		// scoring scheme for aligning profiles
		SWAligner sw_aligner;
		SequenceAlignment local_align( sw_aligner.align( seq1, seq2, css ) );
		std::cout << local_align;
		aligns.push_back(local_align);
		//align_set.insert(local_align);

		NWAligner nw_aligner;
		SequenceAlignment global_align( nw_aligner.align( seq1, seq2, css ) );
		std::cout << global_align;
		aligns.push_back(global_align);
		//align_set.insert(global_align);
		//std::cout << "have " << align_set.size() << " distinct alignments."
		//	<< std::endl;
	}

	std::cout << "done generating alignments." << std::endl;

	//vector1< SequenceAlignment > aligns = align_set.alignments();
	std::cout << "have " << align_set.size() << " distinct alignments."
		<< std::endl;
	for ( vector1< SequenceAlignment >::iterator it = aligns.begin(),
			end = aligns.end();
			it != end; ++it
	) {
		output << *it << std::endl;
	} // for ( it in aligns )

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )

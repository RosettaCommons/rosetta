// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file prof_align.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/L1ScoringScheme.hh>
#include <core/sequence/PairScoringScheme.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/CompositeScoringScheme.hh>

#include <core/scoring/rms_util.hh>

#include <core/pose/Pose.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/comparative_modeling/AlignmentSet.hh>
#include <protocols/jobdist/standard_mains.hh>

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
#include <map>
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


///////////////////////////////////////////////////////////////////////////////

using std::map;
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

int
main( int argc, char* argv [] ) {
	try {

	using core::Real;
	using core::Size;

	// options, random initialization
	devel::init( argc, argv );

	// query and template profiles
	FileName fn1( option[ in::file::pssm ]()[1] );
	FileName fn2( option[ in::file::pssm ]()[2] );

	SequenceProfileOP prof1( new SequenceProfile );
	prof1->read_from_file( fn1 );
	prof1->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()

	SequenceProfileOP prof2( new SequenceProfile );
	prof2->read_from_file( fn2 );
	prof2->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()

	// eliminate leading paths from prof1 and prof2
	prof1->id( FileName( prof1->id() ).base() );
	prof2->id( FileName( prof2->id() ).base() );

	// scoring scheme for aligning profiles
	//std::string const scoring_scheme_type( option[ cm::seq_score ]() );
	std::string const scoring_scheme_type( option[ cm::seq_score ]()[1] );
	ScoringSchemeFactory ssf;
	ScoringSchemeOP ss( ssf.get_scoring_scheme( scoring_scheme_type ) );

	// for the ProfSim scoring scheme, the optimal opening and extension
	// penalties were 2 and 0.2, with a scoring shift of -0.45 applied to
	// all ungapped aligned pairs.
	//Real const max_gap_open    (  option[ cm::max_gap_open ]() );
	Real const min_gap_open    (  option[ cm::min_gap_open ]() );
	//Real const max_gap_extend  (  option[ cm::max_gap_open ]() );
	Real const min_gap_extend  (  option[ cm::min_gap_open ]() );
	//Real const open_step_size  (   0.5 );
	//Real const extend_step_size(   0.1 );

	SWAligner sw_aligner;

	std::string output_fn = option[ out::file::alignment ]();
	utility::io::ozstream output( output_fn );

	AlignmentSet align_set;
	Real const arbitrarily_big_score( 1e5 );
	for ( core::Size i = 1; i <= prof1->length(); ++i ) {
	for ( core::Size j = 1; j <= prof2->length(); ++j ) {
		core::sequence::CompositeScoringSchemeOP css(
			new CompositeScoringScheme()
		);

		css->add_scoring_scheme( ss->clone() );
		css->gap_open  ( min_gap_open   );
		css->gap_extend( min_gap_extend );

		// setup extra-special aligner that guarantees that i,j are aligned
		PairScoringSchemeOP enforcer( new PairScoringScheme() );
		enforcer->add_scored_pair( i, j, arbitrarily_big_score );
		css->add_scoring_scheme( enforcer );

		SequenceAlignment local_align = sw_aligner.align( prof1, prof2, css );
		local_align.score( local_align.score() - arbitrarily_big_score );

		//output << "local_align (" << i << "," << j << ")" << std::endl
		//	<< local_align << std::endl << "--" << std::endl;
		align_set.insert( local_align );
		//std::cout << "done with " << i << "," << j << std::endl;
	} // i
	std::cout << "done with " << i * prof2->length() << " / " << prof1->length() * prof2->length()
		<< ", have " << align_set.size() << " distinct alignments."
		<< std::endl;
	} // j

	using utility::vector1;
	using core::sequence::SequenceAlignment;
	vector1< SequenceAlignment > aligns = align_set.alignments();
	std::cout << "have " << align_set.size() << " distinct alignments." << std::endl;
	for ( vector1< SequenceAlignment >::iterator it = aligns.begin(),
			end = aligns.end();
			it != end; ++it
	) {
		output << *it << std::endl;
	} // for ( it in aligns )

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )

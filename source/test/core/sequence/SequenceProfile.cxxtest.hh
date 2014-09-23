// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/sequence/SequenceProfile.cxxtest.hh
/// @brief
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>

#include <test/core/init_util.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/file/FileName.hh>

// Package Headers

// AUTO-REMOVED #include <core/sequence/util.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>
// AUTO-REMOVED #include <core/sequence/DPScoringScheme.hh>
#include <core/sequence/ProfSimScoringScheme.hh>
#include <core/sequence/NWAligner.hh>

#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/sequence/SWAligner.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

//Auto Headers
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/Aligner.hh>
#include <core/sequence/DP_Matrix.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <boost/functional/hash.hpp>


class SequenceProfile_Tests : public CxxTest::TestSuite {

public:
	SequenceProfile_Tests() {}


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {}

	void test_probs() {
		using core::Real;
		using utility::vector1;
		using namespace core::sequence;

		SequenceProfileOP profile( new SequenceProfile );
		vector1< Real > scores;
		scores.push_back( 1 );
		scores.push_back( 2 );
		scores.push_back( 3 );
		scores.push_back( 4 );
		vector1< vector1< Real > > prof;
		prof.push_back( scores );
		profile->profile( prof );
		profile->convert_profile_to_probs( 1.0 );
		profile->sequence( "X" );
		// Because rescaling creates a new copy internally.
		prof = profile->profile();
		float const TOL(1e-3);
		TS_ASSERT_DELTA( prof[1][4], 0.643, TOL );
		TS_ASSERT_DELTA( prof[1][3], 0.236, TOL );
		TS_ASSERT_DELTA( prof[1][2], 0.087, TOL );
		TS_ASSERT_DELTA( prof[1][1], 0.032, TOL );
	}

	void test_prob_temp() {
		using core::Real;
		using utility::vector1;
		using namespace core::sequence;

		SequenceProfileOP profile( new SequenceProfile );
		vector1< Real > scores;
		scores.push_back( 1 );
		scores.push_back( 2 );
		scores.push_back( 3 );
		scores.push_back( 4 );
		vector1< vector1< Real > > prof;
		prof.push_back( scores );
		profile->profile( prof );
		profile->convert_profile_to_probs( 0.5 );
		profile->sequence( "X" );
		// Because rescaling creates a new copy internally.
		prof = profile->profile();
		float const TOL(1e-3);
		TS_ASSERT_DELTA( prof[1][4], 0.865, TOL );
		TS_ASSERT_DELTA( prof[1][3], 0.117, TOL );
		TS_ASSERT_DELTA( prof[1][2], 0.016, TOL );
		TS_ASSERT_DELTA( prof[1][1], 0.002, TOL );
		TS_ASSERT_DELTA( profile->temp(), 0.5, TOL );
	}

	// ------------------------------------------ //
	/// @brief test reading and writing profile.
	void test_input() {
		using utility::file::FileName;
		using namespace core::sequence;

		SequenceProfileOP prof1( new SequenceProfile( FileName("core/sequence/1aho_.fasta.pssm" ) ) );
		prof1->convert_profile_to_probs(); // was previously implicit in constructor
		SequenceProfileOP prof2( new SequenceProfile( FileName("core/sequence/1aho_hom.fasta.pssm" ) ) );
		prof2->convert_profile_to_probs(); // was previously implicit in constructor
		TS_ASSERT( prof1->length() == 64 );
		TS_ASSERT( prof2->length() == 65 );

//		core::Real const TOLERATED_ERROR( 0.00001 );

		ScoringSchemeOP simple_ss( new SimpleScoringScheme( 6, 1, -8, -1 ) );
		ScoringSchemeOP prof_sim( new ProfSimScoringScheme( -1, -0.1 ) );

		NWAligner nw_align;

		SequenceAlignment simple_align
			= nw_align.align( prof1, prof2, simple_ss );
		SequenceAlignment prof_sim_align
			= nw_align.align( prof1, prof2, prof_sim );

		TS_ASSERT_EQUALS(
			simple_align.sequence(1)->sequence(),
			"VKDGYIVDDVNCTYFCGRNAYCNEECTKLKGESGYCQWASPYGNACYCYKLPDHV--RTKGPGRCH"
		);
		TS_ASSERT_EQUALS(
			simple_align.sequence(2)->sequence(),
			"VRDGYIADDKDCAYFCGRNAYCDEEC-KKGAESGKCWYAGQYGNACWCYKLPDWVPIKQKVSGKCN"
		);

		TS_ASSERT_EQUALS(
			prof_sim_align.sequence(1)->sequence(),
			"VKDGYIVDDVNCTYFCGRNAYCNEECTKLKGESGYCQWASPYGNACYCYKLPDHVRT--KGPGRCH"
		);
		TS_ASSERT_EQUALS(
			prof_sim_align.sequence(2)->sequence(),
			//"VRDGYIADDKDCAYFCGRNAYCDEECK-KGAESGKCWYAGQYGNACWCYKLPDWVPIKQKVSGKCN"
			"VRDGYIADDKDCAYFCGRNAYCDEEC-KKGAESGKCWYAGQYGNACWCYKLPDWVPIKQKVSGKCN"
		);
	} // test_input

	void test_sequence_manipulation() {
		using utility::file::FileName;
		using namespace core::sequence;
		SequenceProfileOP prof( new SequenceProfile( FileName("core/sequence/1aho_.fasta.pssm" ) ) );

		TS_ASSERT( prof->sequence() == "VKDGYIVDDVNCTYFCGRNAYCNEECTKLKGESGYCQWASPYGNACYCYKLPDHVRTKGPGRCH" );
		TS_ASSERT( prof->prof_row(1).size() == 20 );
		prof->delete_position( 1 );
		TS_ASSERT( prof->sequence() == "KDGYIVDDVNCTYFCGRNAYCNEECTKLKGESGYCQWASPYGNACYCYKLPDHVRTKGPGRCH" );
	}

	void test_sequence_entropy() {

	}


	void test_rescale() {
		using core::Real;
		using utility::vector1;
		using namespace core::sequence;

		SequenceProfileOP profile( new SequenceProfile );
		vector1< Real > scores;
		scores.push_back( 1 );
		scores.push_back( 2 );
		scores.push_back( 3 );
		scores.push_back( 4 );
		vector1< vector1< Real > > prof;
		prof.push_back( scores );
		profile->profile( prof );
		profile->sequence( "X" );
		profile->rescale( 0.25 );
		// Because rescale creates a new copy.
		prof = profile->profile();
		float const TOL(1e-3);
		TS_ASSERT_DELTA( prof[1][4], 1.00, TOL );
		TS_ASSERT_DELTA( prof[1][3], 0.75, TOL );
		TS_ASSERT_DELTA( prof[1][2], 0.50, TOL );
		TS_ASSERT_DELTA( prof[1][1], 0.25, TOL );
	}

	void test_global_auto_rescale() {
		using core::Real;
		using utility::vector1;
		using namespace core::sequence;

		SequenceProfileOP profile( new SequenceProfile );
		vector1< Real > scores;
		scores.push_back( 1 );
		scores.push_back( 2 );
		scores.push_back( 3 );
		scores.push_back( 4 );
		vector1< Real > scores2;
		scores2.push_back( -1 );
		scores2.push_back( -2 );
		scores2.push_back( -3 );
		scores2.push_back( 0 );
		vector1< vector1< Real > > prof;
		prof.push_back( scores );
		prof.push_back( scores2 );
		profile->profile( prof );
		profile->sequence( "X" );
		profile->global_auto_rescale();
		// Because rescaliong creates a new copy internally.
		prof = profile->profile();
		float const TOL(1e-3);
		TS_ASSERT_DELTA( prof[1][4], 1.00, TOL );
		TS_ASSERT_DELTA( prof[1][3], 0.75, TOL );
		TS_ASSERT_DELTA( prof[1][2], 0.50, TOL );
		TS_ASSERT_DELTA( prof[1][1], 0.25, TOL );
		TS_ASSERT_DELTA( prof[2][4],  0.00, TOL );
		TS_ASSERT_DELTA( prof[2][3], -0.75, TOL );
		TS_ASSERT_DELTA( prof[2][2], -0.50, TOL );
		TS_ASSERT_DELTA( prof[2][1], -0.25, TOL );
	}

	void test_generate_from_sequence() {
		using core::Real;
		using utility::vector1;
		using namespace core::sequence;

		Sequence seq;
		seq.sequence("TESTSEQ");
		SequenceProfile profile;
		profile.generate_from_sequence( seq, "BLOSUM62" );
		TS_ASSERT( profile.negative_better() == false );
		vector1< vector1< Real > > prof;
		prof = profile.profile();
		float const TOL(1e-3);
		TS_ASSERT_DELTA( prof[1][core::chemical::aa_thr], 5, TOL ); // T->T
		TS_ASSERT_DELTA( prof[1][core::chemical::aa_tyr], -2, TOL ); // T->Y
		TS_ASSERT_DELTA( prof[2][core::chemical::aa_asp], 2, TOL ); // E->D
		TS_ASSERT_DELTA( prof[2][core::chemical::aa_lys], 1, TOL ); // E->K
		TS_ASSERT_DELTA( prof[6][core::chemical::aa_asp], 2, TOL ); // E->D
		TS_ASSERT_DELTA( prof[6][core::chemical::aa_lys], 1, TOL ); // E->K
		TS_ASSERT_DELTA( prof[7][core::chemical::aa_thr], -1, TOL ); // Q->T
		TS_ASSERT_DELTA( prof[7][core::chemical::aa_leu], -2, TOL ); // Q->L
	}
}; // SequenceProfileTests

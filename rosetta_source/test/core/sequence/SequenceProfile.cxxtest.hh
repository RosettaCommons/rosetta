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
// AUTO-REMOVED #include <core/sequence/SWAligner.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

//Auto Headers
#include <core/sequence/SequenceAlignment.hh>
#include <utility/vector1.hh>


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
		profile->scores_to_probs_( scores, 1.0 );
		prof.push_back( scores );
		profile->profile( prof );
		profile->sequence( "X" );
		float const TOL(1e-3);
		TS_ASSERT_DELTA( scores[4], 0.643, TOL );
		TS_ASSERT_DELTA( scores[3], 0.236, TOL );
		TS_ASSERT_DELTA( scores[2], 0.087, TOL );
		TS_ASSERT_DELTA( scores[1], 0.032, TOL );
	}

	// ------------------------------------------ //
	/// @brief test reading and writing profile.
	void test_input() {
		using utility::file::FileName;
		using namespace core::sequence;

		SequenceProfileOP prof1(
			new SequenceProfile( FileName("core/sequence/1aho_.fasta.pssm" ) )
		);
		SequenceProfileOP prof2(
			new SequenceProfile( FileName("core/sequence/1aho_hom.fasta.pssm" ) )
		);
		TS_ASSERT( prof1->length() == 64 );
		TS_ASSERT( prof2->length() == 65 );

//		core::Real const TOLERATED_ERROR( 0.00001 );

		ScoringSchemeOP simple_ss( new SimpleScoringScheme( 6, 1, -8, -1 ) );
		ScoringSchemeOP prof_sim ( new ProfSimScoringScheme( -1, -0.1 ) );

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
		SequenceProfileOP prof(
			new SequenceProfile( FileName("core/sequence/1aho_.fasta.pssm" ) )
		);

		TS_ASSERT( prof->sequence() == "VKDGYIVDDVNCTYFCGRNAYCNEECTKLKGESGYCQWASPYGNACYCYKLPDHVRTKGPGRCH" );
		TS_ASSERT( prof->prof_row(1).size() == 20 );
		prof->delete_position( 1 );
		TS_ASSERT( prof->sequence() == "KDGYIVDDVNCTYFCGRNAYCNEECTKLKGESGYCQWASPYGNACYCYKLPDHVRTKGPGRCH" );
	}

	void test_sequence_entropy() {

	}
}; // SequenceProfileTests

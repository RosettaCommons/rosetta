// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/test_rna_sfxn.cxxtest.hh
/// @brief  Ensure that a pose without full_model_info can be scored as needed
/// by the standard RNA sfxn, just with certain scoreterms zeroed
/// @author Andy Watkins

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

// Project headers
//#include <core/chemical/ResidueTypeSet.hh>
#include <core/types.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>


#include <test/UTracer.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.test_rna_sfxn.cxxtest");

// using declarations
using namespace core;
using namespace scoring;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;

class RNA_ScoreFunctionTest : public CxxTest::TestSuite {

public:
	//chemical::ResidueTypeSetCAP residue_set;

	void setUp() {
		core_init_with_additional_options( "-score:weights rna_res_level_energy4.wts -restore_talaris_behavior" );
	}

	void tearDown() {}

	void test_no_full_model_info_scoring() {
		//one_score_type_test(, "core/scoring/test_in.pdb", "core/scoring/.u");
		Pose pose;
		core::pose::make_pose_from_sequence( pose, "gcgcgcaagcgc", "fa_standard" );
		auto sfxn = get_score_function();

		TS_ASSERT_THROWS_NOTHING( ( *sfxn )( pose ) );

		TS_ASSERT( pose.energies().total_energies()[ other_pose ] == 0 );
		TS_ASSERT( pose.energies().total_energies()[ intermol ] == 0 );
		TS_ASSERT( pose.energies().total_energies()[ loop_close ] == 0 );
		TS_ASSERT( pose.energies().total_energies()[ free_suite ] == 0 );
		TS_ASSERT( pose.energies().total_energies()[ free_2HOprime ] == 0 );
	}

};

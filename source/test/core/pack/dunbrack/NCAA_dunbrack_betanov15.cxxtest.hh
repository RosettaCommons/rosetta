// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/pack/dunbrack/NCAA_dunbrack_betanov15.cxxtest.hh
/// @brief  Test suite for loading custom noncanonical rotamer libraries when using the beta_nov15 scoring function.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

//Movers used for testing:
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/relax/FastRelax.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("core.pack.dunbrack.NCAADunbrackTests_betanov15.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;

class NCAADunbrackTests_betanov15 : public CxxTest::TestSuite {

public:

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options("-beta_nov15");
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	/// @brief Load a pose, mutate one residue to a noncanonical, and
	/// score it (which relies on the fa_dun score, which requires the
	/// rotamer libraries to have been loaded).
	void test_score_with_ncaa()
	{
		TR << "Building pose." << std::endl;
		Pose pose = create_trpcage_ideal_pose();

		TR << "Creating scorefunction." << std::endl;
		ScoreFunctionOP sfxn( core::scoring::get_score_function() );

		TR << "Mutating residue 5 to norvaline." << std::endl;
		//Introduce a norvaline:
		protocols::simple_moves::MutateResidue mutres;
		mutres.set_target(5);
		mutres.set_res_name("NVL");
		mutres.apply(pose);

		TR << "FastRelaxing." << std::endl;
		//Pack and relax
		protocols::relax::FastRelax frlx(sfxn, 1);
		frlx.apply(pose);

		TR << "Final score." << std::endl;
		(*sfxn)(pose);

		TR << "Completed NCAADunbrackTests_betanov15:test_score_with_ncaa()." << std::endl;
	}


};



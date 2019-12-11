// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ArgCationPiEnergy.cxxtest.hh
/// @brief  test suite for ArgCationPiEnergy
/// @author Brian Coventry (bcov@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/types.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>

#include <test/UTracer.hh>
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

using namespace core;


static basic::Tracer TR("core.scoring.methods.ArgCationPiEnergyTests.cxxtest");

class ArgCationPiEnergyTests : public CxxTest::TestSuite
{

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	core::pose::PoseOP
	get_arg_cation_pi_examples() {
		core::pose::PoseOP pose_p = utility::pointer::make_shared<core::pose::Pose>();
		import_pose::pose_from_file( *pose_p, "core/scoring/methods/arg_cation_pi_examples.pdb" );
		return pose_p;
	}

	// We don't actually want to test much here. Just that the values are less than 0
	void test_arg_cation_pi(){
		pose::Pose pose = *get_arg_cation_pi_examples();

		// We only need our scoreterm

		scoring::ScoreFunctionOP sfxn = utility::pointer::make_shared<core::scoring::ScoreFunction>();
		scoring::methods::EnergyMethodOptions options = sfxn->energy_method_options();
		options.arg_cation_pi_his_can_be_pi( true );
		sfxn->set_energy_method_options( options );

		sfxn->set_weight( scoring::arg_cation_pi, 1.0 );

		sfxn->score( pose );

		scoring::Energies const & energies = pose.energies();

		TS_ASSERT( energies.residue_total_energy( 1 ) < 0 ); // trp
		TS_ASSERT( energies.residue_total_energy( 3 ) < 0 ); // tyr
		TS_ASSERT( energies.residue_total_energy( 4 ) < 0 ); // his
		TS_ASSERT( energies.residue_total_energy( 6 ) < 0 ); // phe


		options = sfxn->energy_method_options();
		options.arg_cation_pi_his_can_be_pi( false );
		sfxn->set_energy_method_options( options );

		sfxn->score( pose );

		TS_ASSERT( energies.residue_total_energy( 1 ) < 0 ); // trp
		TS_ASSERT( energies.residue_total_energy( 3 ) < 0 ); // tyr
		TS_ASSERT( energies.residue_total_energy( 4 ) == 0 ); // his
		TS_ASSERT( energies.residue_total_energy( 6 ) < 0 ); // phe


	}

};

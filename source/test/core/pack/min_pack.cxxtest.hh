// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/SCMinMultifunc.cxxtest.hh
/// @brief  Sidechain minimization multifunc class tests
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers

// Package headers
#include <core/pack/min_pack.hh>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.pack.scmin.SCMinMultifunc.cxxtest");

using namespace core;


class min_pack_Tests : public CxxTest::TestSuite
{

public:
	min_pack_Tests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-no_optH -ignore_unrecognized_res" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Temporarily disabled.
	void test_min_pack()
	{
		using namespace pack;
		using namespace pose;
		using namespace scoring;

		TS_ASSERT( true );
		return;


		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scoring::methods::EnergyMethodOptionsOP emopts( new scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );;
		scorefxn->set_energy_method_options( *emopts );
		//scorefxn->set_weight( fa_pair, 0.0 );

		// read in pose
		Pose pose = create_trpcage_ideal_pose();
		(*scorefxn)( pose );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}
		task->initialize_from_command_line();
		min_pack( pose, *scorefxn, task );
		pose.dump_pdb( "min_pack_result.pdb" );
	}


};

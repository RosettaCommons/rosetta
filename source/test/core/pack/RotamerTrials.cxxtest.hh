// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerTrials.cxxtest.hh
/// @brief  test suite for rotamer_trials
/// @author Phil Bradley
/// @author Sergey Lyskov
/// @author P. Douglas Renfrew (renfrew@unc.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include "platform/types.hh"

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/rtmin.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/scoring/ScoreFunction.hh>

#include <core/types.hh>

#include <numeric/angle.functions.hh>

#include <test/UTracer.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.pack.RotamerTrials.cxxtest");

using namespace core;

class RotamerTrials : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;
	Real delta;

public:
	RotamerTrials() {};

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-no_optH" );

		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		// init delta
		delta = 0.01;
	}

	// Shared finalization goes here.
	void tearDown() {
	}


///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------ //
/// @brief test for rotamer_trials
void test_rotamer_trials()
{
	// UTracer log file
	test::UTracer  UT("core/pack/RotamerTrials.u");
	do_rotamer_trials(false, UT);
}

/// @brief test for rotamer trials with minimization
void dont_test_rotamer_trials_with_minimization()
{
	// UTracer log file
	test::UTracer  UT("core/pack/RTMIN.u");
	//do_rotamer_trials(true, UT);
}

void do_rotamer_trials(bool with_minimization, test::UTracer & UT)
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	// init/reset seeds in all RG objects we have to do this inside the test it self function since
	// user could request to run just one singel test.
	core::init::init_random_generators(1101, numeric::random::_RND_TestRun_, "mt19937");

	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );

	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );


	// read in pose
	Pose pose(create_test_in_pdb_pose());
	//core::import_pose::pose_from_pdb( pose, "core/pack/test_in.pdb" );

	// calculate original score
	Energy score_orig = scorefxn( pose );

	// create paker task for rotamer trials
	pack::task::PackerTaskOP task
		( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking();

	// run rotamer trials given our pose, score function and packer task
	if(with_minimization) {
		// Just run RTMIN on a few residues, because it's slow
		for(Size i = 1; i < 100; ++i) task->nonconst_residue_task(i).prevent_repacking();
		for(Size i = 115; i <= pose.total_residue(); ++i) task->nonconst_residue_task(i).prevent_repacking();
		pack::RTMin rtmin;
		rtmin.rtmin( pose, scorefxn, task );
	}
	else pack::rotamer_trials( pose, scorefxn, task );

	// calculate score after rotamer trials
	Energy score_rt = scorefxn( pose );

	// test that scores are lower after running rotamer trials
	TS_ASSERT_LESS_THAN(score_rt, score_orig);

	// test that the score are the same as when the test was created
	// note: last update will sheffler (willsheffler@gmail.com) march 24 2008
	//Energy precomputed_score_orig = 44.2714; //  56.684;
	//Energy precomputed_score_rt = -152.9630;//-155.6311; //-144.5115;
	//TS_ASSERT_DELTA( precomputed_score_orig, score_orig, delta );
	//TS_ASSERT_DELTA( precomputed_score_rt, score_rt, delta );

	UT.abs_tolerance(delta);
	UT << "score_orig = " << score_orig << std::endl;
	UT << "score_rt = " << score_rt << std::endl;

	// test that the rotamers produced are the same as when the test was created
	// note: last update will sheffler (willsheffler@gmail.com) march 24 2008

	Size nres = pose.n_residue();
	for ( Size i = 1; i <= nres; ++i )
		{
			 // get chi angles
			utility::vector1< Real > res_chi = pose.residue( i ).chi();

			// compares values to those in the list
			for ( Size j = 1; j <= res_chi.size(); ++j )
				{
					//Real chi = numeric::nonnegative_principal_angle_degrees( res_chi[ j ] );
					//TR << std::setprecision(3) << std::fixed << std::setw(7) << chi << ", ";
					// std::cerr << "FOOCHI " << chi << std::endl;

					UT << std::setprecision(15);
					//TS_ASSERT_DELTA( precomputed_chi_angles[precompute_counter], chi, delta );

					// fpd
					UT << "residue = " << i << " sin (chi) = " << sin( res_chi[ j ] ) << " cos(chi) = " << cos( res_chi[ j ] ) << std::endl;
				}
		}
};

};

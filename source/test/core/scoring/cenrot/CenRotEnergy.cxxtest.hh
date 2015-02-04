// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/RamachandranEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::RamachandranEnergy.cc
/// @author Andrew Leaver-Fay
/// @author Yuan Liu (wendao@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

// Unit headers
#include <core/scoring/CenRotEnvPairPotential.hh>
#include <core/scoring/methods/CenPairEnergy.hh>
#include <core/scoring/methods/CenRotEnvEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Numeric headers
#include <numeric/conversions.hh>

// Project headers
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

//Auto Headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class CenRotModelEnergyTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		using namespace core;
		core_init_with_additional_options( "-corrections:score:cenrot" );
	}

	// Shared finalization goes here.
	void tearDown() {}

	// env should be tested along with pair
	void test_cen_rot_pair_env_energy()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");
		to_cenrot.apply(pose);
		ScoreFunction sfxn;
		sfxn.set_weight( cen_rot_env, 1.0 );
		sfxn.set_weight( cen_rot_pair, 1.0 );
		sfxn.set_weight( cen_rot_pair_ang, 1.0 );
		Real start_score = sfxn(pose);
		//std::cout.precision(15);
		//std::cout << start_score << std::endl;
		TS_ASSERT_DELTA(start_score, -1.43971912089564, 1e-12);
	}

	void test_cen_rot_vdw()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");
		to_cenrot.apply(pose);
		ScoreFunction sfxn;
		EnergyMethodOptions options(sfxn.energy_method_options());
		options.atom_vdw_atom_type_set_name("centroid_rot:min");
		sfxn.set_energy_method_options(options);
		sfxn.set_weight( vdw, 1.0 );
		Real start_score = sfxn(pose);
		//std::cout.precision(15);
		//std::cout << start_score << std::endl;
		TS_ASSERT_DELTA(start_score, 0.838790208886603, 1e-12);
	}

	void test_cen_rot_dun()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");
		to_cenrot.apply(pose);
		ScoreFunction sfxn;
		sfxn.set_weight( cen_rot_dun, 1.0 );
		Real start_score = sfxn(pose);
		//std::cout.precision(15);
		//std::cout << start_score << std::endl;
		TS_ASSERT_DELTA(start_score, 41.6239745723091, 1e-12);
	}

	void test_cen_rot_repack()
	{
		Pose pose = create_trpcage_ideal_pose();
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");
		to_cenrot.apply(pose);
		ScoreFunctionOP sfop( new ScoreFunction() );
		ScoreFunction & sfxn(*sfop);
		//dun
		sfxn.set_weight( cen_rot_dun, 1.0 );
		//vdw
		EnergyMethodOptions options(sfxn.energy_method_options());
		options.atom_vdw_atom_type_set_name("centroid_rot");
		sfxn.set_energy_method_options(options);
		sfxn.set_weight( vdw, 1.0 );
		//pair and env
		sfxn.set_weight( cen_rot_env, 1.0 );
		sfxn.set_weight( cen_rot_pair, 1.0 );
		sfxn.set_weight( cen_rot_pair_ang, 1.0 );

		//repack
		using namespace core::pack::task;
		TaskFactoryOP main_task_factory( new TaskFactory );
		operation::RestrictToRepackingOP rtrop( new operation::RestrictToRepacking );
		main_task_factory->push_back( rtrop );

		protocols::simple_moves::PackRotamersMover packrotamersmover;
		packrotamersmover.task_factory(main_task_factory);
		packrotamersmover.score_function(sfop);

		packrotamersmover.apply(pose);

		Real final_score = sfxn(pose);
		//std::cout.precision(15);
		//std::cout << final_score << std::endl;
		TS_ASSERT_DELTA( final_score, -11.2961651753465, 1e-12 );
	}

	void test_cen_rot_atomtree_min()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");
		to_cenrot.apply(pose);
		ScoreFunction sfxn;

		//dun
		sfxn.set_weight( cen_rot_dun, 1.0 );
		//vdw
		EnergyMethodOptions options(sfxn.energy_method_options());
		options.atom_vdw_atom_type_set_name("centroid_rot");
		sfxn.set_energy_method_options(options);
		sfxn.set_weight( vdw, 1.0 );
		//pair and env
		sfxn.set_weight( cen_rot_env, 1.0 );
		sfxn.set_weight( cen_rot_pair, 1.0 );
		sfxn.set_weight( cen_rot_pair_ang, 1.0 );

		//atomtree min
		core::kinematics::MoveMap mm;
		mm.set_bb(true); mm.set_chi(true); mm.set_jump(false);
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			core::conformation::Residue const &res_i = pose.residue(ii);
			if ( res_i.aa()!=chemical::aa_gly
			     && res_i.aa()!=chemical::aa_ala
			     && res_i.type().has("CEN")) {
				mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::D ), true );
				mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::THETA ), true );
			}
		}

		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-3, true, false, false );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( 10 );

		core::optimization::AtomTreeMinimizer minimizer;
		minimizer.run( pose, mm, sfxn, minoptions );

		Real start_score = sfxn(pose);
		//std::cout.precision(15);
		//std::cout << start_score << std::endl;
		TS_ASSERT_DELTA(start_score, 36.1618331777519, 1e-5);
	}

	void test_cen_rot_cart_min()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");
		to_cenrot.apply(pose);
		ScoreFunction sfxn;

		//dun
		sfxn.set_weight( cen_rot_dun, 1.0 );
		//vdw
		EnergyMethodOptions options(sfxn.energy_method_options());
		options.atom_vdw_atom_type_set_name("centroid_rot");
		sfxn.set_energy_method_options(options);
		sfxn.set_weight( vdw, 1.0 );
		//pair and env
		sfxn.set_weight( cen_rot_env, 1.0 );
		sfxn.set_weight( cen_rot_pair, 1.0 );
		sfxn.set_weight( cen_rot_pair_ang, 1.0 );
		sfxn.set_weight( cart_bonded, 0.1 );

		//cart min
		core::kinematics::MoveMap mm;
		mm.set_bb(true); mm.set_chi(true); mm.set_jump(false);

		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-3, true, false, false );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( 10 );

		core::optimization::CartesianMinimizer minimizer;
		minimizer.run( pose, mm, sfxn, minoptions );

		Real start_score = sfxn(pose);
		//std::cout.precision(15);
		//std::cout << start_score << std::endl;
		TS_ASSERT_DELTA(start_score, 20.0778485587594, 1e-5);
	}

};



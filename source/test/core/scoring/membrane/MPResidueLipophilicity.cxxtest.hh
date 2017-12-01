// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane.MPResidueLipophilicity.cxxtest.hh
/// @brief  Test the scoring and burial estimations methods of MPResidueLipophilicity
/// @author Jonathan Weinstein

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/membrane/MPResidueLipophilicityEnergy.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

static basic::Tracer TR("core.scoring.membrane.MPResidueLipophilicity.cxxtest");

using namespace core;

class MPResidueLipophilicityTests : public CxxTest::TestSuite {

public:
	MPResidueLipophilicityTests() {};

	void setUp() {
		core_init();
	}

	void tearDown(){}

	void test_burial_sigmoid()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		//core::scoring::ScoreFunction sfxn;
		//sfxn.set_weight( mp_res_lipo, 1 );

		core::Real const TOL(1e-2);
		core::scoring::membrane::MPResidueLipophilicityEnergy mp_res;
		//core::pose::Pose pose = create_1afo_pose();

		// test slope and offset 6
		TS_ASSERT_DELTA( mp_res.burial_sigmoid( 10, 0.15, 20.0 ), 0.8175, TOL);
		TS_ASSERT_DELTA( mp_res.burial_sigmoid( 0, 0.15, 20.0 ), 0.9525, TOL);
		TS_ASSERT_DELTA( mp_res.burial_sigmoid( 5, 0.15, 20.0 ), 0.9046, TOL);
		TS_ASSERT_DELTA( mp_res.burial_sigmoid( 30, 0.15, 20.0 ), 0.1824, TOL);

		// test slope and offset 12
		TS_ASSERT_DELTA( mp_res.burial_sigmoid( 10, 0.5, 475.0 ), 1, TOL);

		// test slope and offset 6 for centorid
		TS_ASSERT_DELTA( mp_res.burial_sigmoid( 10, 0.15, 20.0 ), 0.8175, TOL);

		// test slope and offset 6 for centroid
		TS_ASSERT_DELTA( mp_res.burial_sigmoid( 10, 5, 220.0 ), 1, TOL);
	}

	void test_calc_energy()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose_ );
		//core::pose::Pose pose = create_1afo_pose();

		core::Real const TOL(1e-3);

		core::scoring::membrane::MPResidueLipophilicityEnergy res_lipo;
		EnergyMap emap;
		res_lipo.residue_energy( pose_->residue( 5 ), *pose_, emap );
		TS_ASSERT_DELTA( emap[ MPResidueLipophilicity ], 0.0682, TOL);

		emap[ MPResidueLipophilicity ] = 0;
		res_lipo.residue_energy( pose_->residue( 15 ), *pose_, emap );
		TS_ASSERT_DELTA( emap[ MPResidueLipophilicity ], -3.3328, TOL);

		emap[ MPResidueLipophilicity ] = 0;
		res_lipo.residue_energy( pose_->residue( 20 ), *pose_, emap );
		TS_ASSERT_DELTA( emap[ MPResidueLipophilicity ], 0.2902, TOL);
	}

	// @brief test the centroid level neighbor counting methid centroid_neighbors
	void test_cen_neighbors()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose_ );

		core::Real const TOL(1e0);
		core::pose::PoseOP pose_cen = pose_;
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( *pose_cen );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( MPResidueLipophilicity, 1 );
		sfxn( *pose_ );

		core::scoring::membrane::MPResidueLipophilicityEnergy res_lipo;
		utility::vector1 < core::Size > t = res_lipo.centroid_neighbors( *pose_cen, pose_cen->residue( 5 ) );
		TS_ASSERT_DELTA( res_lipo.centroid_neighbors( *pose_cen, pose_cen->residue( 5 ) )[ 1 ], 29, TOL );
		TS_ASSERT_DELTA( res_lipo.centroid_neighbors( *pose_cen, pose_cen->residue( 10 ) )[ 1 ], 2, TOL );
		TS_ASSERT_DELTA( res_lipo.centroid_neighbors( *pose_cen, pose_cen->residue( 15 ) )[ 1 ], 2, TOL );
	}

	// @brief test teh full atom neighbor counting method neighboring_atoms
	void test_fa_neighbors()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose_ );

		core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
		(*sf)( *pose_ );

		core::Real const TOL(1e-0);

		core::scoring::membrane::MPResidueLipophilicityEnergy res_lipo;
		TS_ASSERT_DELTA( res_lipo.neighboring_atoms( *pose_, pose_->residue( 5 ), 6, 12 )[ 1 ], 44, TOL );
		TS_ASSERT_DELTA( res_lipo.neighboring_atoms( *pose_, pose_->residue( 10 ), 6, 12 )[ 1 ], 12, TOL );
		TS_ASSERT_DELTA( res_lipo.neighboring_atoms( *pose_, pose_->residue( 15 ), 6, 12 )[ 1 ], 18, TOL );
	}

	// @brief test calc_residue_burial
	void test_calc_residue_burial()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose_ );

		core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
		(*sf)( *pose_ );

		core::Real const TOL(1e-3);

		core::scoring::membrane::MPResidueLipophilicityEnergy res_lipo;
		TS_ASSERT_DELTA( res_lipo.calc_residue_burial( *pose_, pose_->residue( 5 ) ), 0.0265, TOL );
		TS_ASSERT_DELTA( res_lipo.calc_residue_burial( *pose_, pose_->residue( 10 ) ), 0.7685, TOL );
		TS_ASSERT_DELTA( res_lipo.calc_residue_burial( *pose_, pose_->residue( 15 ) ), 0.5744, TOL );

		// do the same for a centroid level pose
		core::pose::PoseOP pose_cen = pose_;
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( *pose_cen );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( MPResidueLipophilicity, 1 );
		sfxn( *pose_ );
		TS_ASSERT_DELTA( res_lipo.calc_residue_burial( *pose_cen, pose_cen->residue( 5 ) ), 0.2058, TOL );
		TS_ASSERT_DELTA( res_lipo.calc_residue_burial( *pose_cen, pose_cen->residue( 10 ) ), 0.9370, TOL );
		TS_ASSERT_DELTA( res_lipo.calc_residue_burial( *pose_cen, pose_cen->residue( 15 ) ), 0.9370, TOL );
	}

};

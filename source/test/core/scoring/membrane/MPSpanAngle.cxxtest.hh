// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane.MPSpanAngle.cxxtest.hh
/// @brief  Test the scoring and burial estimations methods of MPSpanAngle
/// @author Jonathan Weinstein

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/membrane/MPSpanAngleEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>

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

static basic::Tracer TR("core.scoring.membrane.MPSpanAngle.cxxtest");

using namespace core;

class MPSpanAngleEnergyTests : public CxxTest::TestSuite {

public:
	MPSpanAngleEnergyTests() {};

	void setUp() {
		core_init();
	}

	void tearDown(){}

	void test_calc_span_score()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::Real const TOL(1e-2);
		core::scoring::membrane::MPSpanAngleEnergy mp_span;

		TS_ASSERT_DELTA( mp_span.calc_ang_score( 0 ), -0.5, TOL);

		TS_ASSERT_DELTA( mp_span.calc_ang_score( 30 ), 1.0299, TOL);

		TS_ASSERT_DELTA( mp_span.calc_ang_score( 60 ), 11.0199, TOL);

		TS_ASSERT_DELTA( mp_span.calc_ang_score( 90 ), 53.7699, TOL);

		TS_ASSERT_DELTA( mp_span.calc_ang_score( 120 ), 153.5799, TOL);
	}

	void test_finalize_total_energy()
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

		core::scoring::membrane::MPSpanAngleEnergy mp_span;
		EnergyMap emap;
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( mp_span_ang, 1 );
		sfxn( *pose_ );

		mp_span.finalize_total_energy( *pose_, sfxn, emap );
		TS_ASSERT_DELTA( emap[ mp_span_ang ], 4.1251, TOL);
	}

	// @brief test the centroid level neighbor counting methid centroid_neighbors
	void test_compute()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose_ );

		core::Real const TOL(1e-3);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( mp_span_ang, 1 );
		sfxn( *pose_ );

		core::scoring::membrane::MPSpanAngleEnergy mp_span;
		utility::vector1< core::Real > res = mp_span.compute( *pose_, TR, true );
		TS_ASSERT_DELTA( res[ 1 ], 0.8139, TOL );
	}

	// @brief test find_helix_vector funciton
	void test_find_helix_vector()
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

		core::scoring::membrane::MPSpanAngleEnergy mp_span;
		TS_ASSERT_DELTA( mp_span.find_helix_vector( *pose_, 33, 58 )[ 1 ][ 1 ], -6.6865, TOL );
		TS_ASSERT_DELTA( mp_span.find_helix_vector( *pose_, 4, 27 )[ 1 ][ 1 ], -3.1370, TOL );
		TS_ASSERT_DELTA( mp_span.find_helix_vector( *pose_, 127, 150 )[ 1 ][ 1 ], 18.2952, TOL );
	}

};

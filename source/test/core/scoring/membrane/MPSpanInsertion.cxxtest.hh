// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane.MPSpanInsertion.cxxtest.hh
/// @brief  Test the scoring and burial estimations methods of MPSpanInsertion
/// @author Jonathan Weinstein

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/membrane/MPSpanInsertionEnergy.hh>
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
#include <core/conformation/membrane/Span.hh>

static basic::Tracer TR("core.scoring.membrane.MPSpanInsertion.cxxtest");

using namespace core;

class MPSpanInsertionEnergyTests : public CxxTest::TestSuite {

public:
	MPSpanInsertionEnergyTests() {};

	void setUp() {
		core_init();
	}

	void tearDown(){}

	// @brief test the calc_span_score method
	void test_calc_span_score()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose_ );

		core::Real const TOL(1e-2);
		core::scoring::membrane::MPSpanInsertionEnergy mp_span;

		TS_ASSERT_DELTA( mp_span.calc_span_score( *pose_, 33, 58 ), 4.2018, TOL);
		TS_ASSERT_DELTA( mp_span.calc_span_score( *pose_, 4, 27 ), 0.5128, TOL);
		TS_ASSERT_DELTA( mp_span.calc_span_score( *pose_, 127, 150 ), -5.6991, TOL);
	}

	// @brief test the spline_by_z method
	void test_spline_by_z()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose_ );

		core::Real const TOL(1e-2);
		core::scoring::membrane::MPSpanInsertionEnergy mp_span;

		TS_ASSERT_DELTA( mp_span.spline_by_z( 'A', 0   ), 0.0000, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'A', 10  ), 0.0000, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'A', 15  ), 0.0000, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'A', -10 ), 0.0000, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'A', -15 ), 0.0000, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'Y', 0   ), 0.8367, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'Y', 10  ), 0.7058, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'Y', 15  ), 0.9544, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'Y', -10 ), 1.4672, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'Y', -15 ), 1.3954, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'W', 0   ), -0.3538, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'W', 10  ), -0.2910, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'W', 15  ), 0.1033, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'W', -10 ), 0.4407, TOL);
		TS_ASSERT_DELTA( mp_span.spline_by_z( 'W', -15 ), 0.7498, TOL);
	}

	void test_finalize_total_energy()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1AFO_AB.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1AFO_AB.span" ) );
		add_memb->apply( *pose_ );

		core::Real const TOL(1e-3);

		core::scoring::membrane::MPSpanInsertionEnergy mp_span;
		EnergyMap emap;
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( mp_span_ang, 1 );
		sfxn( *pose_ );

		mp_span.finalize_total_energy( *pose_, sfxn, emap );
		TS_ASSERT_DELTA( emap[ mp_span_ang ], 0.000, TOL);
	}

	// @brief test the compute method
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

		core::scoring::membrane::MPSpanInsertionEnergy mp_span;
		core::Real res = mp_span.compute( *pose_ );
		TS_ASSERT_DELTA( res, 8.0913, TOL );
	}

	// @brief test cereate updated spans
	void test_create_updated_spans()
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

		core::Real const TOL(1e0);

		core::scoring::membrane::MPSpanInsertionEnergy mp_span;
		utility::vector1< core::conformation::membrane::Span > spans = mp_span.create_updated_span( *pose_ );
		TS_ASSERT_DELTA( spans[ 1 ].start(), 8, TOL );
		TS_ASSERT_DELTA( spans[ 2 ].start(), 37, TOL );
		TS_ASSERT_DELTA( spans[ 3 ].start(), 75, TOL );
		TS_ASSERT_DELTA( spans[ 1 ].end(), 28, TOL );
		TS_ASSERT_DELTA( spans[ 2 ].end(), 57, TOL );
		TS_ASSERT_DELTA( spans[ 3 ].end(), 95, TOL );
		TS_ASSERT_DELTA( spans.size(), 7, TOL );
	}

};

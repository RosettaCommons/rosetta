// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/simple_moves/SixDoFGridDockMoverTests.cxxtest.hh
/// @brief  Unit tests for grid-based sampling in docking
/// @author Odessa Goudy (oda@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/simple_moves/SixDoFGridDockMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <utility/stream_util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("SixDoFGridDockMoverTests");


class SixDoFGridDockMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:
	// typedef protocols::simple_moves::SixDoFGridDockMover SixDoFGridDockMover;
	typedef core::Vector Vector;
	typedef core::Real Real;

public:

	void setUp(){
		core_init();


	}

	void tearDown(){

	}


	void test_xyz_orthogonality() {

		// core::pose::PoseOP pose;
		// pose = core::import_pose::pose_from_file( "core/util/h10_l27_4ayp_grid_a2--3rrq_ha10_3H.pdb");
		// int res1 = 33;
		// int res2a = 122;
		// int res2b = 133;
		//
		//
		// protocols::simple_moves::SixDoFGridDockMover mover;
		// std::tuple< Vector, Vector, Vector > xyz_axes = mover.compute_axes(pose, res1, res2a, res2b);
		//
		// TS_ASSERT_DELTA(axis1.length(), 1, 1e-6);
		// TS_ASSERT_DELTA(axis2.length(), 1, 1e-6);
		// TS_ASSERT_DELTA(axis3.length(), 1, 1e-6);
		//
		// TS_ASSERT_DELTA(axis1.dot(axis2), 0, 1e-6);
		// TS_ASSERT_DELTA(axis2.dot(axis3), 0, 1e-6);
		// TS_ASSERT_DELTA(axis1.dot(axis3), 0, 1e-6);
	}

	void test_parse_range() {
		std::stringstream ss;
		ss << "<SixDoFGridDockMover dof_residue_selector_1=\"resA\" dof_residue_selector_2a=\"resB1\" dof_residue_selector_2b=\"resB2\" name=\"test_range\" range_rot_axis_1=\"-10,10,10\" range_rot_axis_2=\"-10,10,10\" range_rot_axis_3=\"-10,10,10\" range_trans_axis_1=\"-20,20,20\" range_trans_axis_2=\"-20,20,20\" range_trans_axis_3=\"-20,20,20\"/>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		protocols::simple_moves::SixDoFGridDockMover mover;
		utility::vector1 < Real > dof_range = mover.parse_range( tag, "range_rot_axis_1" );
		TS_ASSERT_EQUALS( dof_range[1], -10.0 );
		TS_ASSERT_EQUALS( dof_range[2], 0.0 );
		TS_ASSERT_EQUALS( dof_range[3], 10.0 );
	}

	void test_lex_position() {
		utility::vector1 < utility::vector1 < Real > > dof_values;
		int num = 10;
		// dof_values = { {10, 11, 12}, {20, 21, 22}, {30, 31, ,32}, {40, 41, 42}, {50, 51, 52}, {60, 61, 62} };
		for ( int aa = 1; aa <= 6; ++aa ) {
			utility::vector1< Real > v1;
			for ( int bb = 1; bb <= 2; ++bb ) {
				v1.push_back( num );
				num += 5;
			}
			dof_values.push_back( v1 );
		}
		TR << "dof_values " << dof_values << std::endl;
		protocols::simple_moves::SixDoFGridDockMover mover;
		for ( core::Size ii = 1; ii <= 6; ++ii ) {
			Real value = mover.lex_position( dof_values, 1, ii );
			TR << "dof " << value << std::endl;
			TS_ASSERT_EQUALS( value, dof_values[ ii ][ 1 ] );
		}
		for ( core::Size ii = 1; ii <= 6; ++ii ) {
			Real value = mover.lex_position( dof_values, 64, ii );
			TR << "dof " << value << std::endl;
			TS_ASSERT_EQUALS( value, dof_values[ ii ][ 2 ] );
		}
	}

	void test_pdb_geometries() {
	}

};

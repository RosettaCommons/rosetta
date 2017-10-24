// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file main/source/test/protocols/toolbox/pose_manipulation/pose_manipulation.cxxtest.hh
/// @brief test suite for protocols::toolbox::pose_manipulation::pose_manipulation
/// @author Brian Coventry (bcov@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>

// Package headers
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

// Project headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>

using namespace core;
using namespace protocols::toolbox::pose_manipulation;

static THREAD_LOCAL basic::Tracer TR("protocols.toolbox.pose_manipulation.PoseManipulationTests");


class PoseManipulationTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		trpcage = create_trpcage_ideal_pose();
	}

	void test_repack_functions() {

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		core::pose::Pose pose = trpcage;

		Size target_res = 3;
		char orig_name1 = pose.residue( target_res ).name1();

		core::select::residue_selector::ResidueSubset target_subset( pose.size(), false );
		target_subset[ target_res ] = true;

		repack_this_residue( target_res, pose, scorefxn, true, "" );
		TR << "The residue actually is: " << pose.residue( target_res ).name1() << std::endl;
		TS_ASSERT( orig_name1 == pose.residue( target_res ).name1() );

		repack_these_residues( target_subset, pose, scorefxn, true, "" );
		TR << "The residue actually is: " << pose.residue( target_res ).name1() << std::endl;
		TS_ASSERT( orig_name1 == pose.residue( target_res ).name1() );

		repack_this_residue( target_res, pose, scorefxn, true, "A" );
		TR << "The residue actually is: " << pose.residue( target_res ).name1() << std::endl;
		TS_ASSERT( 'A' == pose.residue( target_res ).name1() );

		repack_these_residues( target_subset, pose, scorefxn, true, "G" );
		TR << "The residue actually is: " << pose.residue( target_res ).name1() << std::endl;
		TS_ASSERT( 'G' == pose.residue( target_res ).name1() );

	}


	void test_linear_rigid_body_moves() {

		core::pose::Pose pose = trpcage;
		core::select::residue_selector::ResidueSubset subset( pose.size(), true );


		numeric::xyzVector<Real> starting_point = pose.residue(1).xyz("N");
		numeric::xyzVector<Real> x_unit( 1, 0, 0 );

		rigid_body_move( x_unit, 0, x_unit, pose, subset );

		numeric::xyzVector<Real> expected_result = starting_point + x_unit;
		TR << "Location is now: " << pose.residue(1).xyz("N").x() << " "
			<< pose.residue(1).xyz("N").y() << " "
			<< pose.residue(1).xyz("N").z() << std::endl;
		TS_ASSERT_DELTA( pose.residue(1).xyz("N").distance( expected_result), 0, 1e-3 );


	}
	void test_rotation_rigid_body_moves() {

		core::pose::Pose pose = trpcage;
		core::select::residue_selector::ResidueSubset subset( pose.size(), true );


		numeric::xyzVector<Real> starting_point = pose.residue(1).xyz("N");
		numeric::xyzVector<Real> x_unit( 1, 0, 0 );
		numeric::xyzVector<Real> y_unit( 0, 1, 0 );
		numeric::xyzVector<Real> z_unit( 0, 0, 1 );
		numeric::xyzVector<Real> origin( 0, 0, 0 );

		rigid_body_move( x_unit, 0, -starting_point+x_unit, pose, subset );

		TR << "Location is now: " << pose.residue(1).xyz("N").x() << " "
			<< pose.residue(1).xyz("N").y() << " "
			<< pose.residue(1).xyz("N").z() << std::endl;

		TS_ASSERT_DELTA( pose.residue(1).xyz("N").distance( x_unit), 0, 1e-3 );

		rigid_body_move( y_unit, -90, origin, pose, subset, origin );

		TR << "Location is now: " << pose.residue(1).xyz("N").x() << " "
			<< pose.residue(1).xyz("N").y() << " "
			<< pose.residue(1).xyz("N").z() << std::endl;

		TS_ASSERT_DELTA( pose.residue(1).xyz("N").distance( z_unit), 0, 1e-3 );


		TR << "If only the next test fails, someone changed pose::center_of_mass()" << std::endl;
		// If this change is correct, you may fix this unit test by setting the correct expected coordinates


		rigid_body_move( x_unit, -90, origin, pose, subset );

		TR << "Location is now: " << pose.residue(1).xyz("N").x() << " "
			<< pose.residue(1).xyz("N").y() << " "
			<< pose.residue(1).xyz("N").z() << std::endl;

		numeric::xyzVector<Real> expected( 0, -14.4915, 6.06065 );
		TS_ASSERT_DELTA( pose.residue(1).xyz("N").distance( expected ), 0, 1e-3 );
	}



private:
	core::pose::Pose trpcage;

};

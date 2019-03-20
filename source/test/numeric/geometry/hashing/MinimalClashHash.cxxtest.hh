// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/geometry/hashing/MinimalClashHash.cxxtest.hh
/// @brief  numeric.geometry.hashing.MinimalClashHash.cxxtest.hh: test suite for MinimalClashHash
/// @author Brian Coventry (bcov@uw.edu)

// Testing headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>


// Project headers
#include <numeric/geometry/hashing/MinimalClashHash.hh>
#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pose/xyzStripeHashPose.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <iostream>
#include <string>
#include <sstream>

//#include <basic/Tracer.hh>

// static basic::Tracer TR( "numeric.geometry.hashing.MinimalClashHash.cxxtest.hh" );

class MinimalClashHashTest : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	void test_that_it_works() {

		core::pose::Pose pose = create_trpcage_ideal_pose();

		core::Size n_heavy_atoms = 0;
		for ( core::Size ires = 1; ires <= pose.size(); ires++ ) {
			n_heavy_atoms += pose.residue(ires).nheavyatoms();
		}

		core::id::AtomID_Map<core::Real> pose_heavy_atoms = core::pose::make_atom_map( pose, core::pose::PoseCoordPickMode_HVY );
		utility::vector1<numeric::geometry::hashing::Ball> pose_balls;
		core::pose::xyzStripeHashPose::extract_pose_balls( pose, pose_balls, pose_heavy_atoms );

		numeric::geometry::hashing::MinimalClashHash clash_hash( 0.25f, 2.0f, pose_balls );

		TS_ASSERT_EQUALS( clash_hash.clash_check_balls( pose_balls ), n_heavy_atoms );

		// Move the pose really far away
		pose.apply_transform_Rx_plus_v( numeric::xyzMatrix< core::Real >::identity(), numeric::xyzVector< core::Real >( 100, 0, 0 ) );

		utility::vector1<numeric::geometry::hashing::Ball> far_balls;
		core::pose::xyzStripeHashPose::extract_pose_balls( pose, far_balls, pose_heavy_atoms );

		TS_ASSERT_EQUALS( clash_hash.clash_check_balls( far_balls ), 0 );
	}


};

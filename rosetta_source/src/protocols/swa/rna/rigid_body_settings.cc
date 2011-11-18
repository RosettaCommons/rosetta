// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNAResidueSampler
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/rigid_body_settings.hh>
#include <protocols/swa/StepWiseUtil.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/tools/make_vector1.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
#include <set>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {

	///////////////////////////////////////////////////////////////
	void
	apply_rigid_body_settings( pose::Pose & pose, utility::vector1< Real > const & rbs, Size const fixed_res, Size const moving_res ){

		using namespace protocols::swa;
		using namespace utility::tools;

		Vector v2, v1, v;
		Matrix M2, M1, M;
		get_base_centroid_and_rotation_matrix( pose, fixed_res, v1, M1 );
		get_base_centroid_and_rotation_matrix( pose, moving_res, v2, M2 );
		create_euler_rotation( M, rbs[1],rbs[2],rbs[3] );
		v = Vector( rbs[4],rbs[5],rbs[6] );

		Matrix M_final = M1 * M * M2.transposed();
		Vector v_final = -1 * ( M1 * ( M * ( M2.transposed() * v2) ) ) + M1 * v + v1;

		rotate( pose, M_final, pose, make_vector1( moving_res ) );
		translate( pose, v_final, pose, make_vector1( moving_res ) );

	}



}
}
}

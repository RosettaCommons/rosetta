// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/util.cc
/// @brief Utility functions for the local_backbone_mover.
/// @author xingjiepan (xingjiepan@gmail.com)

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

// Core headers
#include <core/kinematics/Stub.hh>

// Project headers
#include <protocols/backbone_moves/local_backbone_mover/util.hh>

#include <basic/Tracer.hh>
#include <utility/fixedsizearray1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.util" );


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {


numeric::xyzVector<Real>
xyz_from_internal_coords( numeric::xyzVector<Real> atom1_xyz,
	numeric::xyzVector<Real> atom2_xyz, numeric::xyzVector<Real> atom3_xyz,
	Real phi, Real theta, Real d){
	using numeric::constants::d::pi;

	core::kinematics::Stub stub(atom1_xyz, atom2_xyz, atom3_xyz);

	return stub.v + stub.M * numeric::x_rotation_matrix_radians(phi) \
		* numeric::z_rotation_matrix_radians(pi - theta) \
		* numeric::xyzVector<Real>(d, 0, 0);
}

void
xyz_to_vec1(numeric::xyzVector<Real> const& vec_xyz, fixedsizearray1<Real,3> &vec1){

	for ( Size i=1; i<=3; ++i ) {
		vec1[i] = vec_xyz(i);
	}
}


} //protocols
} //backbone_moves
} //local_backbone_mover



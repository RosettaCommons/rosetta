// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cryst/spacegroup_symmops3.cc
/// @brief Contains helper functions for the symmetry operations.  Split from spacegroups.cc to prevent the
/// 32-bit compilation from running out of memory during compilation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#include <protocols/cryst/spacegroup_symmops3.hh>

//Rosetta core includes:
#include <core/kinematics/RT.hh>

//Rosetta numeric includes:
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>

namespace protocols {
namespace cryst {

void get_symmops_Iminus42m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Iminus42d( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashmmm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashmcc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashnbm__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashnnc__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashmbm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashmnc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashnmm__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashncc__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashmmc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashmcm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashnbc__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashnnm__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashmbc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashmnm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashnmc__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashncm__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_I4slashmmm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(32);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=16; ++ii ) {
		rt_out[16+ii] = rt_out[ii];
		rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_I4slashmcm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(32);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=16; ++ii ) {
		rt_out[16+ii] = rt_out[ii];
		rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_I41slashamd__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(32);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	for ( int ii=1; ii<=16; ++ii ) {
		rt_out[16+ii] = rt_out[ii];
		rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_I41slashacd__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(32);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	for ( int ii=1; ii<=16; ++ii ) {
		rt_out[16+ii] = rt_out[ii];
		rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(3);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
}

void get_symmops_P31( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(3);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
}

void get_symmops_P32( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(3);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
}

void get_symmops_R3__H( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(9);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=3; ++ii ) {
		rt_out[3+ii] = rt_out[ii];
		rt_out[3+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
		rt_out[6+ii] = rt_out[ii];
		rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
}

void get_symmops_Pminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_Rminus3__H( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(18);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=6; ++ii ) {
		rt_out[6+ii] = rt_out[ii];
		rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P312( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P321( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P3112( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P3121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P3212( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P3221( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_R32__H( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(18);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=6; ++ii ) {
		rt_out[6+ii] = rt_out[ii];
		rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P3m1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
}

void get_symmops_P31m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P3c1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
}

void get_symmops_P31c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_R3m__H( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(18);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=6; ++ii ) {
		rt_out[6+ii] = rt_out[ii];
		rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
}

void get_symmops_R3c__H( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(18);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=6; ++ii ) {
		rt_out[6+ii] = rt_out[ii];
		rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
}

void get_symmops_Pminus31m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_Pminus31c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_Pminus3m1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_Pminus3c1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_Rminus3m__H( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(36);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=12; ++ii ) {
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_Rminus3c__H( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(36);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=12; ++ii ) {
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P6( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P61( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P65( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P62( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P64( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P63( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_Pminus6( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(6);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P6slashm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P63slashm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P622( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P6122( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P6522( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P6222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

}
}


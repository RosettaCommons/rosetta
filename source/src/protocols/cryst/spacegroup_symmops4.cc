// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/cryst/spacegroup_symmops4.cc
/// @brief Contains helper functions for the symmetry operations.  Split from spacegroups.cc to prevent the
/// 32-bit compilation from running out of memory during compilation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#include <protocols/cryst/spacegroup_symmops4.hh>

//Rosetta core includes:
#include <core/kinematics/RT.hh>

//Rosetta numeric includes:
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>

namespace protocols {
namespace cryst {

void get_symmops_P6422( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P6322( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
}

void get_symmops_P6mm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P6cc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P63cm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_P63mc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
}

void get_symmops_Pminus6m2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_Pminus6c2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_Pminus62m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_Pminus62c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P6slashmmm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
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
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P6slashmcc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P63slashmcm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P63slashmmc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
}

void get_symmops_P23( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_F23( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=12; ++ii ) {
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[36+ii] = rt_out[ii];
		rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_I23( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=12; ++ii ) {
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_P213( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(12);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_I213( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	for ( int ii=1; ii<=12; ++ii ) {
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Pmminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Pnminus3__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Fmminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(96);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[72+ii] = rt_out[ii];
		rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Fdminus3__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(96);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[72+ii] = rt_out[ii];
		rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Imminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Paminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Iaminus3( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_P432( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_P4232( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_F432( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(96);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[72+ii] = rt_out[ii];
		rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_F4132( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(96);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[72+ii] = rt_out[ii];
		rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_I432( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_P4332( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_P4132( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_I4132( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Pminus43m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Fminus43m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(96);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[72+ii] = rt_out[ii];
		rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Iminus43m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Pminus43n( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(24);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Fminus43c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(96);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[72+ii] = rt_out[ii];
		rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Iminus43d( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	for ( int ii=1; ii<=24; ++ii ) {
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Pmminus3m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Pnminus3n__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Pmminus3n( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Pnminus3m__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(48);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Fmminus3m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(192);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=48; ++ii ) {
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[96+ii] = rt_out[ii];
		rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[144+ii] = rt_out[ii];
		rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Fmminus3c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(192);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	for ( int ii=1; ii<=48; ++ii ) {
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[96+ii] = rt_out[ii];
		rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[144+ii] = rt_out[ii];
		rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Fdminus3m__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(192);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.25) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.75) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.75) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.75) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.25) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.25) );
	for ( int ii=1; ii<=48; ++ii ) {
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[96+ii] = rt_out[ii];
		rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[144+ii] = rt_out[ii];
		rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Fdminus3c__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(192);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.75) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.25) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.75) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.5) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.75) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.75) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.25) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.75) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.25) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.75) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.25) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.5) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.25) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.75) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.25) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.5) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.25) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.25) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.75) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.25) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.75) );
	for ( int ii=1; ii<=48; ++ii ) {
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[96+ii] = rt_out[ii];
		rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[144+ii] = rt_out[ii];
		rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Imminus3m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(96);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=48; ++ii ) {
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_Iaminus3d( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(96);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
	rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
	rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
	for ( int ii=1; ii<=48; ++ii ) {
		rt_out[48+ii] = rt_out[ii];
		rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

void get_symmops_B11m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=2; ++ii ) {
		rt_out[2+ii] = rt_out[ii];
		rt_out[2+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
}

}
}


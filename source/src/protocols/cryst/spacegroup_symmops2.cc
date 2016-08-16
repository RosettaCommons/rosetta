// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cryst/spacegroup_symmops2.cc
/// @brief Contains helper functions for the symmetry operations.  Split from spacegroups.cc to prevent the
/// 32-bit compilation from running out of memory during compilation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#include <protocols/cryst/spacegroup_symmops2.hh>

//Rosetta core includes:
#include <core/kinematics/RT.hh>

//Rosetta numeric includes:
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>

namespace protocols {
namespace cryst {

void get_symmops_Pbca( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pnma( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Cmcm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Cmca( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Cmmm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Cccm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Cmma( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Ccca__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Fmmm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(32);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16+ii] = rt_out[ii];
		rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Fddd__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(32);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16+ii] = rt_out[ii];
		rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24+ii] = rt_out[ii];
		rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Immm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Ibam( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Ibca( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Imma( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P4( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P41( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P42( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P43( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_I4( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_I41( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_Pminus4( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Iminus4( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4slashn__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42slashn__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_I4slashm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_I41slasha__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P422( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4212( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4122( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P41212( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.75) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.25) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P42212( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4322( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P43212( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.75) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.25) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.25) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.75) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_I422( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_I4122( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_P4mm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P4bm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P42cm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P42nm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P4cc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P4nc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P42mc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_P42bc( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_I4mm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_I4cm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_I41md( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_I41cd( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.75) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.25) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
}

void get_symmops_Pminus42m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Pminus42c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Pminus421m( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Pminus421c( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Pminus4m2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Pminus4c2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Pminus4b2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Pminus4n2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Iminus4m2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

void get_symmops_Iminus4c2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=8; ++ii ) {
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
}

}
}


// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cryst/spacegroup_symmops1.cc
/// @brief Contains helper functions for the symmetry operations.  Split from spacegroups.cc to prevent the
/// 32-bit compilation from running out of memory during compilation.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#include <protocols/cryst/spacegroup_symmops1.hh>

//Rosetta core includes:
#include <core/kinematics/RT.hh>

//Rosetta numeric includes:
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>

namespace protocols {
namespace cryst {

void get_symmops_P1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(1);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0, 0 ) );
}

void get_symmops_Pminus1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(2);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(2);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0, 0.5 ) );
}

void get_symmops_P1211( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(2);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0, 0.5 ) );
}

void get_symmops_C121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=2; ++ii ) {
		rt_out[2+ii] = rt_out[ii];
		rt_out[2+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0, 0.5 ) );
}

void get_symmops_P1m1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(2);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0.5, 0 ) );
}

void get_symmops_P1c1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(2);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0.5, 0 ) );
}

void get_symmops_C1m1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=2; ++ii ) {
		rt_out[2+ii] = rt_out[ii];
		rt_out[2+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0.5, 0 ) );
}

void get_symmops_C1c1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=2; ++ii ) {
		rt_out[2+ii] = rt_out[ii];
		rt_out[2+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0.5, 0 ) );
}

void get_symmops_P12slashm1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P121slashm1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_C12slashm1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P12slashc1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P121slashc1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_C12slashc1( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P2221( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P21212( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_P212121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_C2221( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_C222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_F222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_I222( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_I212121( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pmm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pmc21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pcc2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pma2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pca21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pnc2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pmn21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pba2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pna21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pnn2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(4);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Cmm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Cmc21( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Ccc2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Amm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Abm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Ama2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Aba2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Fmm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Fdd2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(16);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8+ii] = rt_out[ii];
		rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12+ii] = rt_out[ii];
		rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Imm2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Iba2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Ima2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	for ( int ii=1; ii<=4; ++ii ) {
		rt_out[4+ii] = rt_out[ii];
		rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	}
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
}

void get_symmops_Pmmm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pnnn__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pccm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pban__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pmma( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pnna( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pmna( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pcca( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pbam( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pccn( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pbcm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pnnm( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pmmn__2( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

void get_symmops_Pbcn( utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc ) {
	rt_out.resize(8);
	rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
	rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
	rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
	rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
	cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
}

}
}


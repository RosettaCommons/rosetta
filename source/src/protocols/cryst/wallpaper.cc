// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <protocols/cryst/wallpaper.hh>
#include <protocols/cryst/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/CrystInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>
#include <cmath>

#include <sstream>
#include <string>
#include <queue>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>



namespace protocols {
namespace cryst {

static basic::Tracer TSW("WallpaperGroup");
///
///
///
WallpaperGroup::WallpaperGroup() {
	name_="";
	setting_=wgNOSETTING;
	a_=b_=0;
	alpha_=90.0;
	V_=0;
	ncopies_=0;
}

WallpaperGroup::WallpaperGroup(std::string name_in) {
	a_=b_=0;
	alpha_=90.0;
	V_=0;
	ncopies_=0;
	set_wallpaper_group( name_in );
}

void
WallpaperGroup::set_wallpaper_group( std::string name_in) {
	// reset lattice
	a_=b_=0;
	alpha_=90.0;
	V_=0;
	ncopies_=0;

	name_ = name_in;
	name_.erase( std::remove( name_.begin(), name_.end(), ' ' ), name_.end() );

	// setting + validation
	if ( name_ == "P1" || name_ == "P112" ) {
		setting_ = wgMONOCLINIC;
	} else if ( name_ == "P121" || name_ == "P2111" || name_ == "C211" || name_ == "P222" || name_ == "P2122" || name_ == "P21212" || name_ == "C222" || name_ == "C211" ) {
		setting_ = wgTETRAGONAL;
	} else if ( name_ == "P4" || name_ == "P422" || name_ == "P4212" ) {
		setting_ = wgCUBIC;
	} else if ( name_ == "P3" || name_ == "P312" || name_ == "P321" || name_ == "P6" || name_ == "P622" ) {
		setting_ = wgHEXAGONAL;
	} else {
		utility_exit_with_message("Unknown WallpaperGroup! "+name_);
	}

	// lookup
	symmops_.clear();
	get_symmops(symmops_, cc_);
}


// sets AND VALIDATES input parameters
void WallpaperGroup::set_parameters(core::Real a_in, core::Real b_in, core::Real alpha_in) {
	if ( setting_ == wgMONOCLINIC ) {
		a_=a_in; b_=b_in; alpha_=alpha_in;
	} else if ( setting_ == wgCUBIC ) {
		a_=a_in; b_=a_in; alpha_=90.0;
	} else if ( setting_ == wgTETRAGONAL ) {
		a_=a_in; b_=b_in; alpha_=90.0;
	} else if ( setting_ == wgHEXAGONAL ) {
		a_=a_in; b_=a_in; alpha_=120.0;
	}

	// transformation matrices
	core::Real ca = cos(DEG2RAD*alpha_);
	core::Real sa = sin(DEG2RAD*alpha_);
	f2c_ = numeric::xyzMatrix<core::Real>::rows(
		a_  , b_ * ca , 0,
		0.0 , b_ * sa , 0,
		0.0 , 0.0     , 1
	);
	c2f_ = numeric::inverse(f2c_);
	V_ = a_*b_* sa;

	// report
	if ( a_!=a_in || b_!=b_in || alpha_!=alpha_in ) {
		TSW << "Overriding input crystal parameters with [ "
			<< a_ << "," << b_ << " , " << alpha_ << " ]" << std::endl;
	}
}


///
///
///
void WallpaperGroup::get_symmops(utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc) const {
	if ( name_ == "P1" ) {
		rt_out.resize(1);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0 ,0 ,0 ) );
	} else if ( name_ == "P112" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0,0,0),numeric::xyzVector<core::Real>(0.5,0.5,0) );
	} else if ( name_ == "P121" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0,0,0),numeric::xyzVector<core::Real>(0.5,0,0.5) );
	} else if ( name_ == "P2111" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5, 0 ) );
	} else if ( name_ == "C211" ) {
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0 ) );
	} else if ( name_ == "P222" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "P2122" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0.5,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "P21212" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "C222" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "P4" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 , 0.5 ,0 ) );
	} else if ( name_ == "P422" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P4212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P3" ) {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	} else if ( name_ == "P312" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0.5 ) );
	} else if ( name_ == "P321" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P6" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0 ) );
	} else if ( name_ == "P622" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	}


	for ( int ii=1; ii<=(int)rt_out.size(); ++ii ) {
		core::Vector T_i = rt_out[ii].get_translation();
		rt_out[ii].set_translation( core::Vector(std::fmod(T_i[0],1), std::fmod(T_i[1],1), std::fmod(T_i[2],1) ) );
	}
}

// name output in pdbheader
std::string WallpaperGroup::pdbname() const {
	if ( name_ == "P1" ) return "P 1";
	if ( name_ == "P112" ) return "P 1 1 2";
	if ( name_ == "P121" ) return "P 1 2 1";
	if ( name_ == "P2111" ) return "P 21 1 1";
	if ( name_ == "C211" ) return "C 2 1 1";
	if ( name_ == "P222" ) return "P 2 2 2";
	if ( name_ == "P2122" ) return "P 21 2 2";
	if ( name_ == "P21212" ) return "P 21 21 2";
	if ( name_ == "C222" ) return "C 2 2 2";
	if ( name_ == "P4" ) return "P 4";
	if ( name_ == "P422" ) return "P 4 2 2";
	if ( name_ == "P4212" ) return "P 4 21 2";
	if ( name_ == "P3" ) return "P 6";
	if ( name_ == "P312" ) return "P 3 1 2";
	if ( name_ == "P321" ) return "P 3 2 1";
	if ( name_ == "P6" ) return "P 6";
	if ( name_ == "P622" ) return "P 6 2 2";
	return "X";
}

}
}


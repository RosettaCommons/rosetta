// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#include <protocols/cryst/util.hh>
#include <protocols/cryst/spacegroup.hh>

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

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
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
#include <math.h>

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


#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323

static THREAD_LOCAL basic::Tracer TSG("spacegroup");

Spacegroup::Spacegroup() {
	name_="";
	setting_=NOSETTING;
	a_=b_=c_=0;
	alpha_=beta_=gamma_=90.0;
	V_=0;
	ncopies_=0;
}

Spacegroup::Spacegroup(std::string name_in) {
	a_=b_=c_=0;
	alpha_=beta_=gamma_=90.0;
	V_=0;
	ncopies_=0;
	set_spacegroup( name_in );
}

void
Spacegroup::set_spacegroup( std::string name_in) {
	// reset lattice
	a_=b_=c_=0;
	alpha_=beta_=gamma_=90.0;
	V_=0;
	ncopies_=0;

	name_ = name_in;
	name_.erase( std::remove( name_.begin(), name_.end(), ' ' ), name_.end() );

	// setting + validation
	if ( name_ == "P1" || name_ == "P-1" ) {
		setting_ = TRICLINIC;
	} else if (
			name_ == "P121" || name_ == "P1211" || name_ == "C121" || name_ == "P1m1" || name_ == "P1c1" || name_ == "C1m1"
			|| name_ == "C1c1" || name_ == "P12/m1" || name_ == "P121/m1" || name_ == "C12/m1" || name_ == "P12/c1"
			|| name_ == "P121/c1" || name_ == "C12/c1"
			) {
		setting_ = MONOCLINIC;
	} else if (
			name_ == "P222" || name_ == "P2221" || name_ == "P21212" || name_ == "P212121" || name_ == "C2221"
			|| name_ == "C222" || name_ == "F222" || name_ == "I222" || name_ == "I212121" || name_ == "Pmm2"
			|| name_ == "Pmc21" || name_ == "Pcc2" || name_ == "Pma2" || name_ == "Pca21" || name_ == "Pnc2"
			|| name_ == "Pmn21" || name_ == "Pba2" || name_ == "Pna21" || name_ == "Pnn2" || name_ == "Cmm2"
			|| name_ == "Cmc21" || name_ == "Ccc2" || name_ == "Amm2" || name_ == "Abm2" || name_ == "Ama2"
			|| name_ == "Aba2" || name_ == "Fmm2" || name_ == "Fdd2" || name_ == "Imm2" || name_ == "Iba2"
			|| name_ == "Ima2" || name_ == "Pmmm" || name_ == "Pnnn:2" || name_ == "Pccm" || name_ == "Pban:2"
			|| name_ == "Pmma" || name_ == "Pnna" || name_ == "Pmna" || name_ == "Pcca" || name_ == "Pbam"
			|| name_ == "Pccn" || name_ == "Pbcm" || name_ == "Pnnm" || name_ == "Pmmn:2" || name_ == "Pbcn"
			|| name_ == "Pbca" || name_ == "Pnma" || name_ == "Cmcm" || name_ == "Cmca" || name_ == "Cmmm"
			|| name_ == "Cccm" || name_ == "Cmma" || name_ == "Ccca:2" || name_ == "Fmmm" || name_ == "Fddd:2"
			|| name_ == "Immm" || name_ == "Ibam" || name_ == "Ibca" || name_ == "Imma"
			) {
		setting_ = ORTHORHOMBIC;
	} else if (
			name_ == "P4" || name_ == "P41" || name_ == "P42" || name_ == "P43" || name_ == "I4" || name_ == "I41"
			|| name_ == "P-4" || name_ == "I-4" || name_ == "P4/m" || name_ == "P42/m" || name_ == "P4/n:2"
			|| name_ == "P42/n:2" || name_ == "I4/m" || name_ == "I41/a:2" || name_ == "P422" || name_ == "P4212"
			|| name_ == "P4122" || name_ == "P41212" || name_ == "P4222" || name_ == "P42212" || name_ == "P4322"
			|| name_ == "P43212" || name_ == "I422" || name_ == "I4122" || name_ == "P4mm" || name_ == "P4bm"
			|| name_ == "P42cm" || name_ == "P42nm" || name_ == "P4cc" || name_ == "P4nc" || name_ == "P42mc"
			|| name_ == "P42bc" || name_ == "I4mm" || name_ == "I4cm" || name_ == "I41md" || name_ == "I41cd"
			|| name_ == "P-42m" || name_ == "P-42c" || name_ == "P-421m" || name_ == "P-421c" || name_ == "P-4m2"
			|| name_ == "P-4c2" || name_ == "P-4b2" || name_ == "P-4n2" || name_ == "I-4m2" || name_ == "I-4c2"
			|| name_ == "I-42m" || name_ == "I-42d" || name_ == "P4/mmm" || name_ == "P4/mcc" || name_ == "P4/nbm:2"
			|| name_ == "P4/nnc:2" || name_ == "P4/mbm" || name_ == "P4/mnc" || name_ == "P4/nmm:2" || name_ == "P4/ncc:2"
			|| name_ == "P42/mmc" || name_ == "P42/mcm" || name_ == "P42/nbc:2" || name_ == "P42/nnm:2" || name_ == "P42/mbc"
			|| name_ == "P42/mnm" || name_ == "P42/nmc:2" || name_ == "P42/ncm:2" || name_ == "I4/mmm" || name_ == "I4/mcm"
			|| name_ == "I41/amd:2" || name_ == "I41/acd:2"
			) {
		setting_ = TETRAGONAL;
	} else if (
			name_ == "P3" || name_ == "P31" || name_ == "P32" || name_ == "R3:H" || name_ == "P-3" || name_ == "R-3:H"
			|| name_ == "P312" || name_ == "P321" || name_ == "P3112" || name_ == "P3121" || name_ == "P3212" || name_ == "P3221"
			|| name_ == "R32:H" || name_ == "P3m1" || name_ == "P31m" || name_ == "P3c1" || name_ == "P31c" || name_ == "R3m:H"
			|| name_ == "R3c:H" || name_ == "P-31m" || name_ == "P-31c" || name_ == "P-3m1" || name_ == "P-3c1"
			|| name_ == "R-3m:H" || name_ == "R-3c:H" || name_ == "P6" || name_ == "P61" || name_ == "P65" || name_ == "P62"
			|| name_ == "P64" || name_ == "P63" || name_ == "P-6" || name_ == "P6/m" || name_ == "P63/m" || name_ == "P622"
			|| name_ == "P6122" || name_ == "P6522" || name_ == "P6222" || name_ == "P6422" || name_ == "P6322"
			|| name_ == "P6mm" || name_ == "P6cc" || name_ == "P63cm" || name_ == "P63mc" || name_ == "P-6m2"
			|| name_ == "P-6c2" || name_ == "P-62m" || name_ == "P-62c" || name_ == "P6/mmm" || name_ == "P6/mcc"
			|| name_ == "P63/mcm" || name_ == "P63/mmc"
			// alt (others?)
			|| name_ == "H3" || name_ == "H32"
			) {
		setting_ = HEXAGONAL;
	} else if (
			name_ == "P23" || name_ == "F23" || name_ == "I23" || name_ == "P213" || name_ == "I213" || name_ == "Pm-3"
			|| name_ == "Pn-3:2" || name_ == "Fm-3" || name_ == "Fd-3:2" || name_ == "Im-3" || name_ == "Pa-3"
			|| name_ == "Ia-3" || name_ == "P432" || name_ == "P4232" || name_ == "F432" || name_ == "F4132"
			|| name_ == "I432" || name_ == "P4332" || name_ == "P4132" || name_ == "I4132" || name_ == "P-43m"
			|| name_ == "F-43m" || name_ == "I-43m" || name_ == "P-43n" || name_ == "F-43c" || name_ == "I-43d"
			|| name_ == "Pm-3m" || name_ == "Pn-3n:2" || name_ == "Pm-3n" || name_ == "Pn-3m:2" || name_ == "Fm-3m"
			|| name_ == "Fm-3c" || name_ == "Fd-3m:2" || name_ == "Fd-3c:2" || name_ == "Im-3m" || name_ == "Ia-3d"
			) {
		setting_ = CUBIC;
	} else {
		utility_exit_with_message("Unknown setting for spacegroup "+name_);
	}

	// lookup
	symmops_.clear();
	get_symmops(symmops_, cc_);
}


// sets AND VALIDATES input parameters
void Spacegroup::set_parameters(core::Real a_in, core::Real b_in, core::Real c_in, core::Real alpha_in, core::Real beta_in, core::Real gamma_in) {
	//TODO: fix this logic for nonstandard settings (e.g., A121)
	if ( setting_ == TRICLINIC ) {
		a_=a_in; b_=b_in; c_=c_in; alpha_=alpha_in; beta_=beta_in; gamma_=gamma_in;
	} else if ( setting_ == MONOCLINIC ) {
		a_=a_in; b_=b_in; c_=c_in; alpha_=90.0; beta_=beta_in; gamma_=90.0;
	} else if ( setting_ == CUBIC ) {
		a_=a_in; b_=a_in; c_=a_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
	} else if ( setting_ == ORTHORHOMBIC ) {
		a_=a_in; b_=b_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
	} else if ( setting_ == TETRAGONAL ) {
		a_=a_in; b_=a_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
	} else if ( setting_ == HEXAGONAL ) {
		a_=a_in; b_=a_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=120.0;
	}

	// transformation matrices
	core::Real ca = cos(DEG2RAD*alpha_), cb = cos(DEG2RAD*beta_), cg = cos(DEG2RAD*gamma_);
	core::Real /*sa = sin(DEG2RAD*alpha_),*/ sb = sin(DEG2RAD*beta_), sg = sin(DEG2RAD*gamma_);
	f2c_ = numeric::xyzMatrix<core::Real>::rows(
		a_  , b_ * cg , c_ * cb,
		0.0 , b_ * sg , c_ * (ca - cb*cg) / sg,
		0.0 , 0.0     , c_ * sb * sqrt(1.0 - numeric::square((cb*cg - ca)/(sb*sg)))
	);
	c2f_ = numeric::inverse(f2c_);
	V_ = a_*b_*c_* sqrt(1-numeric::square(ca)-numeric::square(cb)-numeric::square(cg)+2*ca*cb*cg);

	// report
	if ( a_!=a_in || b_!=b_in || c_!=c_in || alpha_!=alpha_in || beta_!=beta_in || gamma_!=gamma_in ) {
		TSG << "Overriding input crystal parameters ["
			<< a_in << "," << b_in << "," << c_in << " , " << alpha_in << ","  << beta_in << ","  << gamma_in << "] "
			<< "] with [ " << a_ << "," << b_ << "," << c_ << " , " << alpha_ << ","  << beta_ << ","  << gamma_ << " ]" << std::endl;
	}
}


///
///
///
void Spacegroup::get_symmops(utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc) const {
	if ( name_ == "P1" ) {
		rt_out.resize(1);;
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0 ,0 ,0 ) );
	}
	if ( name_ == "P-1" ) {
		rt_out.resize(2);;
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P121" || name_ == "P2" ) {
		rt_out.resize(2);;
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0,0,0),numeric::xyzVector<core::Real>(0.5,0,0.5) );
	}
	if ( name_ == "P1211" || name_ == "P21" ) {
		rt_out.resize(2);;
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0 ,0.5 ) );
	}
	if ( name_ == "C121" || name_ == "C2" ) {
		rt_out.resize(4);;
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0 ,0.5 ) );
	}
	if ( name_ == "P1m1" ) {
		rt_out.resize(306);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=18; ++ii ) {
			rt_out[18+ii] = rt_out[ii];
			rt_out[18+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[54+ii] = rt_out[ii];
			rt_out[54+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[90+ii] = rt_out[ii];
			rt_out[90+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[108+ii] = rt_out[ii];
			rt_out[108+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[126+ii] = rt_out[ii];
			rt_out[126+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[162+ii] = rt_out[ii];
			rt_out[162+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[180+ii] = rt_out[ii];
			rt_out[180+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[198+ii] = rt_out[ii];
			rt_out[198+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[234+ii] = rt_out[ii];
			rt_out[234+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[252+ii] = rt_out[ii];
			rt_out[252+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[270+ii] = rt_out[ii];
			rt_out[270+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0.5, 0 ) );
	}
	if ( name_ == "P1c1" ) {
		rt_out.resize(18);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=6; ++ii ) {
			rt_out[6+ii] = rt_out[ii];
			rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0.5, 0 ) );
	}
	if ( name_ == "C1m1" ) {
		rt_out.resize(180);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=18; ++ii ) {
			rt_out[18+ii] = rt_out[ii];
			rt_out[18+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[54+ii] = rt_out[ii];
			rt_out[54+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[90+ii] = rt_out[ii];
			rt_out[90+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[108+ii] = rt_out[ii];
			rt_out[108+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[126+ii] = rt_out[ii];
			rt_out[126+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[162+ii] = rt_out[ii];
			rt_out[162+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0.5, 0 ) );
	}
	if ( name_ == "C1c1" ) {
		rt_out.resize(324);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=18; ++ii ) {
			rt_out[18+ii] = rt_out[ii];
			rt_out[18+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[54+ii] = rt_out[ii];
			rt_out[54+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[90+ii] = rt_out[ii];
			rt_out[90+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[108+ii] = rt_out[ii];
			rt_out[108+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[126+ii] = rt_out[ii];
			rt_out[126+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[162+ii] = rt_out[ii];
			rt_out[162+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[180+ii] = rt_out[ii];
			rt_out[180+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[198+ii] = rt_out[ii];
			rt_out[198+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[234+ii] = rt_out[ii];
			rt_out[234+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[252+ii] = rt_out[ii];
			rt_out[252+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[270+ii] = rt_out[ii];
			rt_out[270+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[306+ii] = rt_out[ii];
			rt_out[306+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0, 0.5, 0 ) );
	}
	if ( name_ == "P12/m1" ) {
		rt_out.resize(1330);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=38; ++ii ) {
			rt_out[38+ii] = rt_out[ii];
			rt_out[38+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[76+ii] = rt_out[ii];
			rt_out[76+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[114+ii] = rt_out[ii];
			rt_out[114+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[152+ii] = rt_out[ii];
			rt_out[152+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[190+ii] = rt_out[ii];
			rt_out[190+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[228+ii] = rt_out[ii];
			rt_out[228+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[266+ii] = rt_out[ii];
			rt_out[266+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[304+ii] = rt_out[ii];
			rt_out[304+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[342+ii] = rt_out[ii];
			rt_out[342+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[380+ii] = rt_out[ii];
			rt_out[380+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[418+ii] = rt_out[ii];
			rt_out[418+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[456+ii] = rt_out[ii];
			rt_out[456+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[494+ii] = rt_out[ii];
			rt_out[494+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[532+ii] = rt_out[ii];
			rt_out[532+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[570+ii] = rt_out[ii];
			rt_out[570+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[608+ii] = rt_out[ii];
			rt_out[608+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[646+ii] = rt_out[ii];
			rt_out[646+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[684+ii] = rt_out[ii];
			rt_out[684+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[722+ii] = rt_out[ii];
			rt_out[722+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[760+ii] = rt_out[ii];
			rt_out[760+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[798+ii] = rt_out[ii];
			rt_out[798+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[836+ii] = rt_out[ii];
			rt_out[836+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[874+ii] = rt_out[ii];
			rt_out[874+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[912+ii] = rt_out[ii];
			rt_out[912+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[950+ii] = rt_out[ii];
			rt_out[950+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[988+ii] = rt_out[ii];
			rt_out[988+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1026+ii] = rt_out[ii];
			rt_out[1026+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[1064+ii] = rt_out[ii];
			rt_out[1064+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1102+ii] = rt_out[ii];
			rt_out[1102+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[1140+ii] = rt_out[ii];
			rt_out[1140+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1178+ii] = rt_out[ii];
			rt_out[1178+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[1216+ii] = rt_out[ii];
			rt_out[1216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1254+ii] = rt_out[ii];
			rt_out[1254+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[1292+ii] = rt_out[ii];
			rt_out[1292+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P121/m1" ) {
		rt_out.resize(36);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "C12/m1" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P12/c1" ) {
		rt_out.resize(612);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=36; ++ii ) {
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[108+ii] = rt_out[ii];
			rt_out[108+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[180+ii] = rt_out[ii];
			rt_out[180+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[252+ii] = rt_out[ii];
			rt_out[252+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[324+ii] = rt_out[ii];
			rt_out[324+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[360+ii] = rt_out[ii];
			rt_out[360+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[396+ii] = rt_out[ii];
			rt_out[396+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[432+ii] = rt_out[ii];
			rt_out[432+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[468+ii] = rt_out[ii];
			rt_out[468+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[504+ii] = rt_out[ii];
			rt_out[504+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[540+ii] = rt_out[ii];
			rt_out[540+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[576+ii] = rt_out[ii];
			rt_out[576+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P121/c1" ) {
		rt_out.resize(324);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		for ( int ii=1; ii<=36; ++ii ) {
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[108+ii] = rt_out[ii];
			rt_out[108+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[180+ii] = rt_out[ii];
			rt_out[180+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[252+ii] = rt_out[ii];
			rt_out[252+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "C12/c1" ) {
		rt_out.resize(360);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=36; ++ii ) {
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[108+ii] = rt_out[ii];
			rt_out[108+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[180+ii] = rt_out[ii];
			rt_out[180+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[252+ii] = rt_out[ii];
			rt_out[252+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[324+ii] = rt_out[ii];
			rt_out[324+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P222" ) {
		rt_out.resize(2520);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[49] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[50] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[51] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[52] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[53] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[54] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[55] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[56] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[57] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[58] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[59] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[60] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[61] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[62] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[63] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[64] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[65] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[66] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[67] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[68] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[69] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[70] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[71] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[72] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=72; ++ii ) {
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[360+ii] = rt_out[ii];
			rt_out[360+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[432+ii] = rt_out[ii];
			rt_out[432+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[504+ii] = rt_out[ii];
			rt_out[504+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[576+ii] = rt_out[ii];
			rt_out[576+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[648+ii] = rt_out[ii];
			rt_out[648+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[720+ii] = rt_out[ii];
			rt_out[720+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[792+ii] = rt_out[ii];
			rt_out[792+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[864+ii] = rt_out[ii];
			rt_out[864+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[936+ii] = rt_out[ii];
			rt_out[936+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[1008+ii] = rt_out[ii];
			rt_out[1008+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1080+ii] = rt_out[ii];
			rt_out[1080+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[1152+ii] = rt_out[ii];
			rt_out[1152+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1224+ii] = rt_out[ii];
			rt_out[1224+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[1296+ii] = rt_out[ii];
			rt_out[1296+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1368+ii] = rt_out[ii];
			rt_out[1368+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[1440+ii] = rt_out[ii];
			rt_out[1440+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1512+ii] = rt_out[ii];
			rt_out[1512+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[1584+ii] = rt_out[ii];
			rt_out[1584+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1656+ii] = rt_out[ii];
			rt_out[1656+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[1728+ii] = rt_out[ii];
			rt_out[1728+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1800+ii] = rt_out[ii];
			rt_out[1800+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[1872+ii] = rt_out[ii];
			rt_out[1872+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1944+ii] = rt_out[ii];
			rt_out[1944+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[2016+ii] = rt_out[ii];
			rt_out[2016+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[2088+ii] = rt_out[ii];
			rt_out[2088+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[2160+ii] = rt_out[ii];
			rt_out[2160+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[2232+ii] = rt_out[ii];
			rt_out[2232+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[2304+ii] = rt_out[ii];
			rt_out[2304+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[2376+ii] = rt_out[ii];
			rt_out[2376+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[2448+ii] = rt_out[ii];
			rt_out[2448+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P2221" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P21212" ) {
		rt_out.resize(36);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P212121" ) {
		rt_out.resize(36);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "C2221" ) {
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
	if ( name_ == "C222" ) {
		rt_out.resize(72);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[60+ii] = rt_out[ii];
			rt_out[60+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "F222" ) {
		rt_out.resize(96);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[60+ii] = rt_out[ii];
			rt_out[60+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[84+ii] = rt_out[ii];
			rt_out[84+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "I222" ) {
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
	if ( name_ == "I212121" ) {
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
	if ( name_ == "Pmm2" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pmc21" ) {
		rt_out.resize(36);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pcc2" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pma2" ) {
		rt_out.resize(36);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pca21" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pnc2" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pmn21" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pba2" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pna21" ) {
		rt_out.resize(36);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pnn2" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Cmm2" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Cmc21" ) {
		rt_out.resize(72);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[60+ii] = rt_out[ii];
			rt_out[60+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Ccc2" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[168+ii] = rt_out[ii];
			rt_out[168+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[264+ii] = rt_out[ii];
			rt_out[264+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Amm2" ) {
		rt_out.resize(72);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[60+ii] = rt_out[ii];
			rt_out[60+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Abm2" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[168+ii] = rt_out[ii];
			rt_out[168+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[264+ii] = rt_out[ii];
			rt_out[264+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Ama2" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[168+ii] = rt_out[ii];
			rt_out[168+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[264+ii] = rt_out[ii];
			rt_out[264+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Aba2" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[168+ii] = rt_out[ii];
			rt_out[168+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[264+ii] = rt_out[ii];
			rt_out[264+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Fmm2" ) {
		rt_out.resize(336);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[168+ii] = rt_out[ii];
			rt_out[168+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[216+ii] = rt_out[ii];
			rt_out[216+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[264+ii] = rt_out[ii];
			rt_out[264+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[312+ii] = rt_out[ii];
			rt_out[312+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Fdd2" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[60+ii] = rt_out[ii];
			rt_out[60+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[84+ii] = rt_out[ii];
			rt_out[84+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[108+ii] = rt_out[ii];
			rt_out[108+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[132+ii] = rt_out[ii];
			rt_out[132+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Imm2" ) {
		rt_out.resize(120);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[60+ii] = rt_out[ii];
			rt_out[60+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[84+ii] = rt_out[ii];
			rt_out[84+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[108+ii] = rt_out[ii];
			rt_out[108+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Iba2" ) {
		rt_out.resize(72);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[60+ii] = rt_out[ii];
			rt_out[60+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Ima2" ) {
		rt_out.resize(72);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[60+ii] = rt_out[ii];
			rt_out[60+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0 ) );
	}
	if ( name_ == "Pmmm" ) {
		rt_out.resize(308);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=28; ++ii ) {
			rt_out[28+ii] = rt_out[ii];
			rt_out[28+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[56+ii] = rt_out[ii];
			rt_out[56+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[84+ii] = rt_out[ii];
			rt_out[84+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[112+ii] = rt_out[ii];
			rt_out[112+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[140+ii] = rt_out[ii];
			rt_out[140+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[168+ii] = rt_out[ii];
			rt_out[168+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[196+ii] = rt_out[ii];
			rt_out[196+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[224+ii] = rt_out[ii];
			rt_out[224+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[252+ii] = rt_out[ii];
			rt_out[252+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[280+ii] = rt_out[ii];
			rt_out[280+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pnnn:2" ) {
		rt_out.resize(32);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		for ( int ii=1; ii<=16; ++ii ) {
			rt_out[16+ii] = rt_out[ii];
			rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pccm" ) {
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
	if ( name_ == "Pban:2" ) {
		rt_out.resize(128);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[64+ii] = rt_out[ii];
			rt_out[64+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pmma" ) {
		rt_out.resize(200);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=40; ++ii ) {
			rt_out[40+ii] = rt_out[ii];
			rt_out[40+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[80+ii] = rt_out[ii];
			rt_out[80+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[160+ii] = rt_out[ii];
			rt_out[160+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pnna" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pmna" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pcca" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pbam" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pccn" ) {
		rt_out.resize(72);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pbcm" ) {
		rt_out.resize(72);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pnnm" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pmmn:2" ) {
		rt_out.resize(128);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[64+ii] = rt_out[ii];
			rt_out[64+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pbcn" ) {
		rt_out.resize(200);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=40; ++ii ) {
			rt_out[40+ii] = rt_out[ii];
			rt_out[40+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[80+ii] = rt_out[ii];
			rt_out[80+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[160+ii] = rt_out[ii];
			rt_out[160+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pbca" ) {
		rt_out.resize(288);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Pnma" ) {
		rt_out.resize(32);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		for ( int ii=1; ii<=16; ++ii ) {
			rt_out[16+ii] = rt_out[ii];
			rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Cmcm" ) {
		rt_out.resize(336);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Cmca" ) {
		rt_out.resize(576);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[336+ii] = rt_out[ii];
			rt_out[336+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[384+ii] = rt_out[ii];
			rt_out[384+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[432+ii] = rt_out[ii];
			rt_out[432+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[480+ii] = rt_out[ii];
			rt_out[480+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[528+ii] = rt_out[ii];
			rt_out[528+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Cmmm" ) {
		rt_out.resize(576);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[336+ii] = rt_out[ii];
			rt_out[336+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[384+ii] = rt_out[ii];
			rt_out[384+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[432+ii] = rt_out[ii];
			rt_out[432+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[480+ii] = rt_out[ii];
			rt_out[480+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[528+ii] = rt_out[ii];
			rt_out[528+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Cccm" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Cmma" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Ccca:2" ) {
		rt_out.resize(784);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[49] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[50] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[51] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[52] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[53] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[54] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[55] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[56] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=56; ++ii ) {
			rt_out[56+ii] = rt_out[ii];
			rt_out[56+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[112+ii] = rt_out[ii];
			rt_out[112+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[168+ii] = rt_out[ii];
			rt_out[168+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[224+ii] = rt_out[ii];
			rt_out[224+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[280+ii] = rt_out[ii];
			rt_out[280+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[336+ii] = rt_out[ii];
			rt_out[336+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[392+ii] = rt_out[ii];
			rt_out[392+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[448+ii] = rt_out[ii];
			rt_out[448+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[504+ii] = rt_out[ii];
			rt_out[504+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[560+ii] = rt_out[ii];
			rt_out[560+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[616+ii] = rt_out[ii];
			rt_out[616+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[672+ii] = rt_out[ii];
			rt_out[672+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[728+ii] = rt_out[ii];
			rt_out[728+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Fmmm" ) {
		rt_out.resize(2112);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[49] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[50] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[51] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[52] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[53] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[54] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[55] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[56] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[57] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[58] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[59] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[60] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[61] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[62] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[63] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[64] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[65] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[66] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[67] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[68] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[69] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[70] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[71] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[72] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[73] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[74] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[75] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[76] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[77] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[78] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[79] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[80] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[81] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[82] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[83] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[84] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[85] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[86] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[87] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[88] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=88; ++ii ) {
			rt_out[88+ii] = rt_out[ii];
			rt_out[88+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[176+ii] = rt_out[ii];
			rt_out[176+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[264+ii] = rt_out[ii];
			rt_out[264+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[352+ii] = rt_out[ii];
			rt_out[352+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[440+ii] = rt_out[ii];
			rt_out[440+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[528+ii] = rt_out[ii];
			rt_out[528+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[616+ii] = rt_out[ii];
			rt_out[616+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[704+ii] = rt_out[ii];
			rt_out[704+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[792+ii] = rt_out[ii];
			rt_out[792+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[880+ii] = rt_out[ii];
			rt_out[880+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[968+ii] = rt_out[ii];
			rt_out[968+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[1056+ii] = rt_out[ii];
			rt_out[1056+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1144+ii] = rt_out[ii];
			rt_out[1144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[1232+ii] = rt_out[ii];
			rt_out[1232+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1320+ii] = rt_out[ii];
			rt_out[1320+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[1408+ii] = rt_out[ii];
			rt_out[1408+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1496+ii] = rt_out[ii];
			rt_out[1496+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[1584+ii] = rt_out[ii];
			rt_out[1584+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1672+ii] = rt_out[ii];
			rt_out[1672+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[1760+ii] = rt_out[ii];
			rt_out[1760+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[1848+ii] = rt_out[ii];
			rt_out[1848+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[1936+ii] = rt_out[ii];
			rt_out[1936+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[2024+ii] = rt_out[ii];
			rt_out[2024+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Fddd:2" ) {
		rt_out.resize(128);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		for ( int ii=1; ii<=16; ++ii ) {
			rt_out[16+ii] = rt_out[ii];
			rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[64+ii] = rt_out[ii];
			rt_out[64+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[80+ii] = rt_out[ii];
			rt_out[80+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[112+ii] = rt_out[ii];
			rt_out[112+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Immm" ) {
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
	if ( name_ == "Ibam" ) {
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
	if ( name_ == "Ibca" ) {
		rt_out.resize(144);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[120+ii] = rt_out[ii];
			rt_out[120+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Imma" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		for ( int ii=1; ii<=16; ++ii ) {
			rt_out[16+ii] = rt_out[ii];
			rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "P4" ) {
		rt_out.resize(484);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=44; ++ii ) {
			rt_out[44+ii] = rt_out[ii];
			rt_out[44+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[88+ii] = rt_out[ii];
			rt_out[88+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[132+ii] = rt_out[ii];
			rt_out[132+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[176+ii] = rt_out[ii];
			rt_out[176+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[220+ii] = rt_out[ii];
			rt_out[220+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[264+ii] = rt_out[ii];
			rt_out[264+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[308+ii] = rt_out[ii];
			rt_out[308+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[352+ii] = rt_out[ii];
			rt_out[352+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[396+ii] = rt_out[ii];
			rt_out[396+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[440+ii] = rt_out[ii];
			rt_out[440+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
	}
	if ( name_ == "P41" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
	}
	if ( name_ == "P42" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
	}
	if ( name_ == "P43" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0 ) );
	}
	if ( name_ == "I4" ) {
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
	if ( name_ == "I41" ) {
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
	if ( name_ == "P-4" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "I-4" ) {
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
	if ( name_ == "P4/m" ) {
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
	if ( name_ == "P42/m" ) {
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
	if ( name_ == "P4/n:2" ) {
		rt_out.resize(32);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		for ( int ii=1; ii<=16; ++ii ) {
			rt_out[16+ii] = rt_out[ii];
			rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P42/n:2" ) {
		rt_out.resize(32);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		for ( int ii=1; ii<=16; ++ii ) {
			rt_out[16+ii] = rt_out[ii];
			rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "I4/m" ) {
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
	if ( name_ == "I41/a:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		for ( int ii=1; ii<=16; ++ii ) {
			rt_out[16+ii] = rt_out[ii];
			rt_out[16+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P422" ) {
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
	if ( name_ == "P4212" ) {
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
	if ( name_ == "P4122" ) {
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
	if ( name_ == "P41212" ) {
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
	if ( name_ == "P4222" ) {
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
	if ( name_ == "P42212" ) {
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
	if ( name_ == "P4322" ) {
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
	if ( name_ == "P43212" ) {
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
	if ( name_ == "I422" ) {
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
	if ( name_ == "I4122" ) {
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
	if ( name_ == "P4mm" ) {
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
	if ( name_ == "P4bm" ) {
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
	if ( name_ == "P42cm" ) {
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
	if ( name_ == "P42nm" ) {
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
	if ( name_ == "P4cc" ) {
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
	if ( name_ == "P4nc" ) {
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
	if ( name_ == "P42mc" ) {
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
	if ( name_ == "P42bc" ) {
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
	if ( name_ == "I4mm" ) {
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
	if ( name_ == "I4cm" ) {
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
	if ( name_ == "I41md" ) {
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
	if ( name_ == "I41cd" ) {
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
	if ( name_ == "P-42m" ) {
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
	if ( name_ == "P-42c" ) {
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
	if ( name_ == "P-421m" ) {
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
	if ( name_ == "P-421c" ) {
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
	if ( name_ == "P-4m2" ) {
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
	if ( name_ == "P-4c2" ) {
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
	if ( name_ == "P-4b2" ) {
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
	if ( name_ == "P-4n2" ) {
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
	if ( name_ == "I-4m2" ) {
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
	if ( name_ == "I-4c2" ) {
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
	if ( name_ == "I-42m" ) {
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
	if ( name_ == "I-42d" ) {
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
	if ( name_ == "P4/mmm" ) {
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
	if ( name_ == "P4/mcc" ) {
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
	if ( name_ == "P4/nbm:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P4/nnc:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P4/mbm" ) {
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
	if ( name_ == "P4/mnc" ) {
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
	if ( name_ == "P4/nmm:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P4/ncc:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P42/mmc" ) {
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
	if ( name_ == "P42/mcm" ) {
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
	if ( name_ == "P42/nbc:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P42/nnm:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P42/mbc" ) {
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
	if ( name_ == "P42/mnm" ) {
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
	if ( name_ == "P42/nmc:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P42/ncm:2" ) {
		rt_out.resize(64);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "I4/mmm" ) {
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
	if ( name_ == "I4/mcm" ) {
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
	if ( name_ == "I41/amd:2" ) {
		rt_out.resize(128);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[64+ii] = rt_out[ii];
			rt_out[64+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "I41/acd:2" ) {
		rt_out.resize(128);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.75) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.25) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.25) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.75) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		for ( int ii=1; ii<=32; ++ii ) {
			rt_out[32+ii] = rt_out[ii];
			rt_out[32+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
			rt_out[64+ii] = rt_out[ii];
			rt_out[64+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 0.5, 0.5 ) );
	}
	if ( name_ == "P3" ) {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
	}
	if ( name_ == "P31" ) {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
	}
	if ( name_ == "P32" ) {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
	}
	if ( name_ == "R3:H" || name_ == "H3" ) {
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
	if ( name_ == "P-3" ) {
		rt_out.resize(18);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=9; ++ii ) {
			rt_out[9+ii] = rt_out[ii];
			rt_out[9+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
	}
	if ( name_ == "R-3:H" ) {
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
	if ( name_ == "P312" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
	}
	if ( name_ == "P321" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
	}
	if ( name_ == "P3112" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
	}
	if ( name_ == "P3121" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
	}
	if ( name_ == "P3212" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
	}
	if ( name_ == "P3221" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
	}
	if ( name_ == "R32:H" || name_ == "H32" ) {
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
	if ( name_ == "P3m1" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
	}
	if ( name_ == "P31m" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
	}
	if ( name_ == "P3c1" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
	}
	if ( name_ == "P31c" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
	}
	if ( name_ == "R3m:H" ) {
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
	if ( name_ == "R3c:H" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
			rt_out[36+ii] = rt_out[ii];
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0 ) );
	}
	if ( name_ == "P-31m" ) {
		rt_out.resize(36);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=18; ++ii ) {
			rt_out[18+ii] = rt_out[ii];
			rt_out[18+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
	}
	if ( name_ == "P-31c" ) {
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
	if ( name_ == "P-3m1" ) {
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
	if ( name_ == "P-3c1" ) {
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
	if ( name_ == "R-3m:H" ) {
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
	if ( name_ == "R-3c:H" ) {
		rt_out.resize(96);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
			rt_out[72+ii] = rt_out[ii];
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0.5 ) );
	}
	if ( name_ == "P6" ) {
		rt_out.resize(36);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=18; ++ii ) {
			rt_out[18+ii] = rt_out[ii];
			rt_out[18+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
	}
	if ( name_ == "P61" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
	}
	if ( name_ == "P65" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
	}
	if ( name_ == "P62" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
	}
	if ( name_ == "P64" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
	}
	if ( name_ == "P63" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 0 ) );
	}
	if ( name_ == "P-6" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.666666666666667, 0.666666666666667, 0.5 ) );
	}
	if ( name_ == "P6/m" ) {
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
	if ( name_ == "P63/m" ) {
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
	if ( name_ == "P622" ) {
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
	if ( name_ == "P6122" ) {
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
	if ( name_ == "P6522" ) {
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
	if ( name_ == "P6222" ) {
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
	if ( name_ == "P6422" ) {
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
	if ( name_ == "P6322" ) {
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
	if ( name_ == "P6mm" ) {
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
	if ( name_ == "P6cc" ) {
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
	if ( name_ == "P63cm" ) {
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
	if ( name_ == "P63mc" ) {
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
	if ( name_ == "P-6m2" ) {
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
	if ( name_ == "P-6c2" ) {
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
	if ( name_ == "P-62m" ) {
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
	if ( name_ == "P-62c" ) {
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
	if ( name_ == "P6/mmm" ) {
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
	if ( name_ == "P6/mcc" ) {
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
	if ( name_ == "P63/mcm" ) {
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
	if ( name_ == "P63/mmc" ) {
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
	if ( name_ == "P23" ) {
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
	if ( name_ == "F23" ) {
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
	if ( name_ == "I23" ) {
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
	if ( name_ == "P213" ) {
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
	if ( name_ == "I213" ) {
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
	if ( name_ == "Pm-3" ) {
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
	if ( name_ == "Pn-3:2" ) {
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
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
	}
	if ( name_ == "Fm-3" ) {
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
	if ( name_ == "Fd-3:2" ) {
		rt_out.resize(384);
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
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
		for ( int ii=1; ii<=48; ++ii ) {
			rt_out[48+ii] = rt_out[ii];
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[144+ii] = rt_out[ii];
			rt_out[144+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[240+ii] = rt_out[ii];
			rt_out[240+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[336+ii] = rt_out[ii];
			rt_out[336+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Im-3" ) {
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
	if ( name_ == "Pa-3" ) {
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
	if ( name_ == "Ia-3" ) {
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
	if ( name_ == "P432" ) {
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
	if ( name_ == "P4232" ) {
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
	if ( name_ == "F432" ) {
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
	if ( name_ == "F4132" ) {
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
	if ( name_ == "I432" ) {
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
	if ( name_ == "P4332" ) {
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
	if ( name_ == "P4132" ) {
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
	if ( name_ == "I4132" ) {
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
	if ( name_ == "P-43m" ) {
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
	if ( name_ == "F-43m" ) {
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
	if ( name_ == "I-43m" ) {
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
	if ( name_ == "P-43n" ) {
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
	if ( name_ == "F-43c" ) {
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
	if ( name_ == "I-43d" ) {
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
	if ( name_ == "Pm-3m" ) {
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
	if ( name_ == "Pn-3n:2" ) {
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
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[49] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[50] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[51] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[52] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[53] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[54] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[55] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[56] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[57] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[58] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[59] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[60] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[61] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[62] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[63] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[64] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[65] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[66] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[67] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[68] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[69] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[70] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[71] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[72] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[73] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[74] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[75] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[76] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[77] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[78] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[79] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[80] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[81] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[82] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[83] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[84] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[85] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[86] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[87] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[88] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[89] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[90] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[91] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[92] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[93] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[94] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[95] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[96] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		for ( int ii=1; ii<=96; ++ii ) {
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
	}
	if ( name_ == "Pm-3n" ) {
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
	if ( name_ == "Pn-3m:2" ) {
		rt_out.resize(192);
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
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[49] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[50] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[51] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[52] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[53] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[54] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[55] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[56] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[57] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[58] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[59] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[60] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[61] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[62] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[63] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[64] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[65] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[66] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[67] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[68] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[69] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[70] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[71] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[72] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[73] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[74] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[75] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[76] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[77] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[78] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[79] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[80] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[81] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[82] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[83] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[84] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[85] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[86] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[87] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[88] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[89] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[90] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[91] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[92] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[93] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[94] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[95] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[96] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		for ( int ii=1; ii<=96; ++ii ) {
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(1, 1, 1 ) );
	}
	if ( name_ == "Fm-3m" ) {
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
	if ( name_ == "Fm-3c" ) {
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
	if ( name_ == "Fd-3m:2" ) {
		rt_out.resize(768);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[49] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[50] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[51] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
		rt_out[52] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		rt_out[53] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[54] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
		rt_out[55] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		rt_out[56] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[57] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[58] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[59] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
		rt_out[60] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		rt_out[61] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[62] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
		rt_out[63] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		rt_out[64] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[65] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[66] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.25) );
		rt_out[67] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
		rt_out[68] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.75) );
		rt_out[69] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[70] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
		rt_out[71] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[72] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.75) );
		rt_out[73] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[74] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[75] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
		rt_out[76] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[77] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[78] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
		rt_out[79] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[80] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[81] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[82] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[83] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
		rt_out[84] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[85] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[86] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
		rt_out[87] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.25) );
		rt_out[88] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[89] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[90] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.75) );
		rt_out[91] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
		rt_out[92] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.25) );
		rt_out[93] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[94] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
		rt_out[95] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.75) );
		rt_out[96] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.25) );
		for ( int ii=1; ii<=96; ++ii ) {
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[384+ii] = rt_out[ii];
			rt_out[384+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[480+ii] = rt_out[ii];
			rt_out[480+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[576+ii] = rt_out[ii];
			rt_out[576+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[672+ii] = rt_out[ii];
			rt_out[672+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Fd-3c:2" ) {
		rt_out.resize(768);
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
		rt_out[25] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[26] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[27] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[28] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[29] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[30] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[31] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[32] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[33] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[34] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[35] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[36] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[37] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[38] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[39] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[40] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[41] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[42] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		rt_out[43] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[44] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[45] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[46] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[47] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[48] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[49] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[50] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.25) );
		rt_out[51] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
		rt_out[52] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.75) );
		rt_out[53] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[54] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.5) );
		rt_out[55] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.75) );
		rt_out[56] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[57] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[58] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.25) );
		rt_out[59] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0.5) );
		rt_out[60] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.75) );
		rt_out[61] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.25) );
		rt_out[62] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0.5) );
		rt_out[63] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.75) );
		rt_out[64] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[65] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[66] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.25,0.75) );
		rt_out[67] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.75,0) );
		rt_out[68] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.25) );
		rt_out[69] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[70] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.25,0) );
		rt_out[71] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0,0.75) );
		rt_out[72] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.25) );
		rt_out[73] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,-1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[74] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, -1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.75) );
		rt_out[75] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,1,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
		rt_out[76] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 1,0,0, 0,0,-1 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.25) );
		rt_out[77] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[78] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, -1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.5) );
		rt_out[79] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,-1,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.25) );
		rt_out[80] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 1,0,0, 0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[81] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, -1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[82] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,-1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.75) );
		rt_out[83] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 1,0,0, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0.5) );
		rt_out[84] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,1, 0,-1,0 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.25) );
		rt_out[85] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.75,0.75) );
		rt_out[86] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( -1,0,0, 0,0,-1, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0.5) );
		rt_out[87] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, -1,0,0, 0,1,0 )  , numeric::xyzVector<core::Real>(0.75,0.5,0.25) );
		rt_out[88] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 1,0,0, 0,0,1, 0,1,0 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[89] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,-1, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[90] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,-1,0, 0,0,1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.75,0.25) );
		rt_out[91] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,-1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.25,0) );
		rt_out[92] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,-1, 1,0,0 )  , numeric::xyzVector<core::Real>(0.25,0.5,0.75) );
		rt_out[93] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,1,0, 1,0,0 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[94] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,1,0, 0,0,1, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0.75,0) );
		rt_out[95] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,-1, 0,1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0.75,0,0.25) );
		rt_out[96] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows( 0,0,1, 0,-1,0, -1,0,0 )  , numeric::xyzVector<core::Real>(0,0.25,0.75) );
		for ( int ii=1; ii<=96; ++ii ) {
			rt_out[96+ii] = rt_out[ii];
			rt_out[96+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[192+ii] = rt_out[ii];
			rt_out[192+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[288+ii] = rt_out[ii];
			rt_out[288+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
			rt_out[384+ii] = rt_out[ii];
			rt_out[384+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0,0) );
			rt_out[480+ii] = rt_out[ii];
			rt_out[480+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[576+ii] = rt_out[ii];
			rt_out[576+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[672+ii] = rt_out[ii];
			rt_out[672+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0), numeric::xyzVector<core::Real>(0.5, 0.5, 0.5 ) );
	}
	if ( name_ == "Im-3m" ) {
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
	if ( name_ == "Ia-3d" ) {
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

	// cenops can make translation outside of [0,1], correct this
	for ( int ii=1; ii<=(int)rt_out.size(); ++ii ) {
		core::Vector T_i = rt_out[ii].get_translation();
		rt_out[ii].set_translation( core::Vector(std::fmod(T_i[0],1), std::fmod(T_i[1],1), std::fmod(T_i[2],1) ) );
	}
}

// name output in pdbheader
std::string Spacegroup::pdbname() const {
	if ( name_ == "P1" ) return "P 1";
	if ( name_ == "P-1" ) return "P -1";
	if ( name_ == "P121" ) return "P 1 2 1";
	if ( name_ == "P1211" ) return "P 1 21 1";
	if ( name_ == "C121" ) return "C 1 2 1";
	if ( name_ == "P1m1" ) return "P 1 m 1";
	if ( name_ == "P1c1" ) return "P 1 c 1";
	if ( name_ == "C1m1" ) return "C 1 m 1";
	if ( name_ == "C1c1" ) return "C 1 c 1";
	if ( name_ == "P12/m1" ) return "P 1 2/m 1";
	if ( name_ == "P121/m1" ) return "P 1 21/m 1";
	if ( name_ == "C12/m1" ) return "C 1 2/m 1";
	if ( name_ == "P12/c1" ) return "P 1 2/c 1";
	if ( name_ == "P121/c1" ) return "P 1 21/c 1";
	if ( name_ == "C12/c1" ) return "C 1 2/c 1";
	if ( name_ == "P222" ) return "P 2 2 2";
	if ( name_ == "P2221" ) return "P 2 2 21";
	if ( name_ == "P21212" ) return "P 21 21 2";
	if ( name_ == "P212121" ) return "P 21 21 21";
	if ( name_ == "C2221" ) return "C 2 2 21";
	if ( name_ == "C222" ) return "C 2 2 2";
	if ( name_ == "F222" ) return "F 2 2 2";
	if ( name_ == "I222" ) return "I 2 2 2";
	if ( name_ == "I212121" ) return "I 21 21 21";
	if ( name_ == "Pmm2" ) return "P m m 2";
	if ( name_ == "Pmc21" ) return "P m c 21";
	if ( name_ == "Pcc2" ) return "P c c 2";
	if ( name_ == "Pma2" ) return "P m a 2";
	if ( name_ == "Pca21" ) return "P c a 21";
	if ( name_ == "Pnc2" ) return "P n c 2";
	if ( name_ == "Pmn21" ) return "P m n 21";
	if ( name_ == "Pba2" ) return "P b a 2";
	if ( name_ == "Pna21" ) return "P n a 21";
	if ( name_ == "Pnn2" ) return "P n n 2";
	if ( name_ == "Cmm2" ) return "C m m 2";
	if ( name_ == "Cmc21" ) return "C m c 21";
	if ( name_ == "Ccc2" ) return "C c c 2";
	if ( name_ == "Amm2" ) return "A m m 2";
	if ( name_ == "Abm2" ) return "A b m 2";
	if ( name_ == "Ama2" ) return "A m a 2";
	if ( name_ == "Aba2" ) return "A b a 2";
	if ( name_ == "Fmm2" ) return "F m m 2";
	if ( name_ == "Fdd2" ) return "F d d 2";
	if ( name_ == "Imm2" ) return "I m m 2";
	if ( name_ == "Iba2" ) return "I b a 2";
	if ( name_ == "Ima2" ) return "I m a 2";
	if ( name_ == "Pmmm" ) return "P m m m";
	if ( name_ == "Pnnn:2" ) return "P n n n :2";
	if ( name_ == "Pccm" ) return "P c c m";
	if ( name_ == "Pban:2" ) return "P b a n :2";
	if ( name_ == "Pmma" ) return "P m m a";
	if ( name_ == "Pnna" ) return "P n n a";
	if ( name_ == "Pmna" ) return "P m n a";
	if ( name_ == "Pcca" ) return "P c c a";
	if ( name_ == "Pbam" ) return "P b a m";
	if ( name_ == "Pccn" ) return "P c c n";
	if ( name_ == "Pbcm" ) return "P b c m";
	if ( name_ == "Pnnm" ) return "P n n m";
	if ( name_ == "Pmmn:2" ) return "P m m n :2";
	if ( name_ == "Pbcn" ) return "P b c n";
	if ( name_ == "Pbca" ) return "P b c a";
	if ( name_ == "Pnma" ) return "P n m a";
	if ( name_ == "Cmcm" ) return "C m c m";
	if ( name_ == "Cmca" ) return "C m c a";
	if ( name_ == "Cmmm" ) return "C m m m";
	if ( name_ == "Cccm" ) return "C c c m";
	if ( name_ == "Cmma" ) return "C m m a";
	if ( name_ == "Ccca:2" ) return "C c c a :2";
	if ( name_ == "Fmmm" ) return "F m m m";
	if ( name_ == "Fddd:2" ) return "F d d d :2";
	if ( name_ == "Immm" ) return "I m m m";
	if ( name_ == "Ibam" ) return "I b a m";
	if ( name_ == "Ibca" ) return "I b c a";
	if ( name_ == "Imma" ) return "I m m a";
	if ( name_ == "P4" ) return "P 4";
	if ( name_ == "P41" ) return "P 41";
	if ( name_ == "P42" ) return "P 42";
	if ( name_ == "P43" ) return "P 43";
	if ( name_ == "I4" ) return "I 4";
	if ( name_ == "I41" ) return "I 41";
	if ( name_ == "P-4" ) return "P -4";
	if ( name_ == "I-4" ) return "I -4";
	if ( name_ == "P4/m" ) return "P 4/m";
	if ( name_ == "P42/m" ) return "P 42/m";
	if ( name_ == "P4/n:2" ) return "P 4/n :2";
	if ( name_ == "P42/n:2" ) return "P 42/n :2";
	if ( name_ == "I4/m" ) return "I 4/m";
	if ( name_ == "I41/a:2" ) return "I 41/a :2";
	if ( name_ == "P422" ) return "P 4 2 2";
	if ( name_ == "P4212" ) return "P 4 21 2";
	if ( name_ == "P4122" ) return "P 41 2 2";
	if ( name_ == "P41212" ) return "P 41 21 2";
	if ( name_ == "P4222" ) return "P 42 2 2";
	if ( name_ == "P42212" ) return "P 42 21 2";
	if ( name_ == "P4322" ) return "P 43 2 2";
	if ( name_ == "P43212" ) return "P 43 21 2";
	if ( name_ == "I422" ) return "I 4 2 2";
	if ( name_ == "I4122" ) return "I 41 2 2";
	if ( name_ == "P4mm" ) return "P 4 m m";
	if ( name_ == "P4bm" ) return "P 4 b m";
	if ( name_ == "P42cm" ) return "P 42 c m";
	if ( name_ == "P42nm" ) return "P 42 n m";
	if ( name_ == "P4cc" ) return "P 4 c c";
	if ( name_ == "P4nc" ) return "P 4 n c";
	if ( name_ == "P42mc" ) return "P 42 m c";
	if ( name_ == "P42bc" ) return "P 42 b c";
	if ( name_ == "I4mm" ) return "I 4 m m";
	if ( name_ == "I4cm" ) return "I 4 c m";
	if ( name_ == "I41md" ) return "I 41 m d";
	if ( name_ == "I41cd" ) return "I 41 c d";
	if ( name_ == "P-42m" ) return "P -4 2 m";
	if ( name_ == "P-42c" ) return "P -4 2 c";
	if ( name_ == "P-421m" ) return "P -4 21 m";
	if ( name_ == "P-421c" ) return "P -4 21 c";
	if ( name_ == "P-4m2" ) return "P -4 m 2";
	if ( name_ == "P-4c2" ) return "P -4 c 2";
	if ( name_ == "P-4b2" ) return "P -4 b 2";
	if ( name_ == "P-4n2" ) return "P -4 n 2";
	if ( name_ == "I-4m2" ) return "I -4 m 2";
	if ( name_ == "I-4c2" ) return "I -4 c 2";
	if ( name_ == "I-42m" ) return "I -4 2 m";
	if ( name_ == "I-42d" ) return "I -4 2 d";
	if ( name_ == "P4/mmm" ) return "P 4/m m m";
	if ( name_ == "P4/mcc" ) return "P 4/m c c";
	if ( name_ == "P4/nbm:2" ) return "P 4/n b m :2";
	if ( name_ == "P4/nnc:2" ) return "P 4/n n c :2";
	if ( name_ == "P4/mbm" ) return "P 4/m b m";
	if ( name_ == "P4/mnc" ) return "P 4/m n c";
	if ( name_ == "P4/nmm:2" ) return "P 4/n m m :2";
	if ( name_ == "P4/ncc:2" ) return "P 4/n c c :2";
	if ( name_ == "P42/mmc" ) return "P 42/m m c";
	if ( name_ == "P42/mcm" ) return "P 42/m c m";
	if ( name_ == "P42/nbc:2" ) return "P 42/n b c :2";
	if ( name_ == "P42/nnm:2" ) return "P 42/n n m :2";
	if ( name_ == "P42/mbc" ) return "P 42/m b c";
	if ( name_ == "P42/mnm" ) return "P 42/m n m";
	if ( name_ == "P42/nmc:2" ) return "P 42/n m c :2";
	if ( name_ == "P42/ncm:2" ) return "P 42/n c m :2";
	if ( name_ == "I4/mmm" ) return "I 4/m m m";
	if ( name_ == "I4/mcm" ) return "I 4/m c m";
	if ( name_ == "I41/amd:2" ) return "I 41/a m d :2";
	if ( name_ == "I41/acd:2" ) return "I 41/a c d :2";
	if ( name_ == "P3" ) return "P 3";
	if ( name_ == "P31" ) return "P 31";
	if ( name_ == "P32" ) return "P 32";
	if ( name_ == "R3:H" || name_ == "H3" ) return "R 3 :H";
	if ( name_ == "P-3" ) return "P -3";
	if ( name_ == "R-3:H" ) return "R -3 :H";
	if ( name_ == "P312" ) return "P 3 1 2";
	if ( name_ == "P321" ) return "P 3 2 1";
	if ( name_ == "P3112" ) return "P 31 1 2";
	if ( name_ == "P3121" ) return "P 31 2 1";
	if ( name_ == "P3212" ) return "P 32 1 2";
	if ( name_ == "P3221" ) return "P 32 2 1";
	if ( name_ == "R32:H" || name_ == "H32" ) return "R 3 2 :H";
	if ( name_ == "P3m1" ) return "P 3 m 1";
	if ( name_ == "P31m" ) return "P 3 1 m";
	if ( name_ == "P3c1" ) return "P 3 c 1";
	if ( name_ == "P31c" ) return "P 3 1 c";
	if ( name_ == "R3m:H" ) return "R 3 m :H";
	if ( name_ == "R3c:H" ) return "R 3 c :H";
	if ( name_ == "P-31m" ) return "P -3 1 m";
	if ( name_ == "P-31c" ) return "P -3 1 c";
	if ( name_ == "P-3m1" ) return "P -3 m 1";
	if ( name_ == "P-3c1" ) return "P -3 c 1";
	if ( name_ == "R-3m:H" ) return "R -3 m :H";
	if ( name_ == "R-3c:H" ) return "R -3 c :H";
	if ( name_ == "P6" ) return "P 6";
	if ( name_ == "P61" ) return "P 61";
	if ( name_ == "P65" ) return "P 65";
	if ( name_ == "P62" ) return "P 62";
	if ( name_ == "P64" ) return "P 64";
	if ( name_ == "P63" ) return "P 63";
	if ( name_ == "P-6" ) return "P -6";
	if ( name_ == "P6/m" ) return "P 6/m";
	if ( name_ == "P63/m" ) return "P 63/m";
	if ( name_ == "P622" ) return "P 6 2 2";
	if ( name_ == "P6122" ) return "P 61 2 2";
	if ( name_ == "P6522" ) return "P 65 2 2";
	if ( name_ == "P6222" ) return "P 62 2 2";
	if ( name_ == "P6422" ) return "P 64 2 2";
	if ( name_ == "P6322" ) return "P 63 2 2";
	if ( name_ == "P6mm" ) return "P 6 m m";
	if ( name_ == "P6cc" ) return "P 6 c c";
	if ( name_ == "P63cm" ) return "P 63 c m";
	if ( name_ == "P63mc" ) return "P 63 m c";
	if ( name_ == "P-6m2" ) return "P -6 m 2";
	if ( name_ == "P-6c2" ) return "P -6 c 2";
	if ( name_ == "P-62m" ) return "P -6 2 m";
	if ( name_ == "P-62c" ) return "P -6 2 c";
	if ( name_ == "P6/mmm" ) return "P 6/m m m";
	if ( name_ == "P6/mcc" ) return "P 6/m c c";
	if ( name_ == "P63/mcm" ) return "P 63/m c m";
	if ( name_ == "P63/mmc" ) return "P 63/m m c";
	if ( name_ == "P23" ) return "P 2 3";
	if ( name_ == "F23" ) return "F 2 3";
	if ( name_ == "I23" ) return "I 2 3";
	if ( name_ == "P213" ) return "P 21 3";
	if ( name_ == "I213" ) return "I 21 3";
	if ( name_ == "Pm-3" ) return "P m -3";
	if ( name_ == "Pn-3:2" ) return "P n -3 :2";
	if ( name_ == "Fm-3" ) return "F m -3";
	if ( name_ == "Fd-3:2" ) return "F d -3 :2";
	if ( name_ == "Im-3" ) return "I m -3";
	if ( name_ == "Pa-3" ) return "P a -3";
	if ( name_ == "Ia-3" ) return "I a -3";
	if ( name_ == "P432" ) return "P 4 3 2";
	if ( name_ == "P4232" ) return "P 42 3 2";
	if ( name_ == "F432" ) return "F 4 3 2";
	if ( name_ == "F4132" ) return "F 41 3 2";
	if ( name_ == "I432" ) return "I 4 3 2";
	if ( name_ == "P4332" ) return "P 43 3 2";
	if ( name_ == "P4132" ) return "P 41 3 2";
	if ( name_ == "I4132" ) return "I 41 3 2";
	if ( name_ == "P-43m" ) return "P -4 3 m";
	if ( name_ == "F-43m" ) return "F -4 3 m";
	if ( name_ == "I-43m" ) return "I -4 3 m";
	if ( name_ == "P-43n" ) return "P -4 3 n";
	if ( name_ == "F-43c" ) return "F -4 3 c";
	if ( name_ == "I-43d" ) return "I -4 3 d";
	if ( name_ == "Pm-3m" ) return "P m -3 m";
	if ( name_ == "Pn-3n:2" ) return "P n -3 n :2";
	if ( name_ == "Pm-3n" ) return "P m -3 n";
	if ( name_ == "Pn-3m:2" ) return "P n -3 m :2";
	if ( name_ == "Fm-3m" ) return "F m -3 m";
	if ( name_ == "Fm-3c" ) return "F m -3 c";
	if ( name_ == "Fd-3m:2" ) return "F d -3 m :2";
	if ( name_ == "Fd-3c:2" ) return "F d -3 c :2";
	if ( name_ == "Im-3m" ) return "I m -3 m";
	if ( name_ == "Ia-3d" ) return "I a -3 d";

	return "X";
}

}
}

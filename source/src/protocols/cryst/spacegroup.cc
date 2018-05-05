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

#include <protocols/cryst/util.hh>
#include <protocols/cryst/spacegroup.hh>
#include <protocols/cryst/spacegroup_symmops1.hh>
#include <protocols/cryst/spacegroup_symmops2.hh>
#include <protocols/cryst/spacegroup_symmops3.hh>
#include <protocols/cryst/spacegroup_symmops4.hh>

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


#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323

static basic::Tracer TSG("spacegroup");

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
			|| name_ == "P121/c1" || name_ == "P121/n1" || name_ == "231" || name_ == "C12/c1"
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
	//if ( a_!=a_in || b_!=b_in || c_!=c_in || alpha_!=alpha_in || beta_!=beta_in || gamma_!=gamma_in ) {
	TSG << "Overriding input crystal parameters [ "
		<< a_in << "," << b_in << "," << c_in << " , " << alpha_in << ","  << beta_in << ","  << gamma_in << " ] "
		<< "with [ " << a_ << "," << b_ << "," << c_ << " , " << alpha_ << ","  << beta_ << ","  << gamma_ << " ]" << std::endl;
	//}
}


///
///
///
void Spacegroup::get_symmops(utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc) const {
	if ( name_ == "P1" ) { get_symmops_P1(rt_out, cc ); }
	if ( name_ == "P-1" ) { get_symmops_Pminus1(rt_out, cc ); }
	if ( name_ == "P121" || name_ == "P2" ) { get_symmops_P121(rt_out, cc ); }
	if ( name_ == "P1211" || name_ == "P21" ) { get_symmops_P1211(rt_out, cc ); }
	if ( name_ == "C121" ) { get_symmops_C121(rt_out, cc ); }
	if ( name_ == "P1m1" ) { get_symmops_P1m1(rt_out, cc ); }
	if ( name_ == "P1c1" ) { get_symmops_P1c1(rt_out, cc ); }
	if ( name_ == "C1m1" ) { get_symmops_C1m1(rt_out, cc ); }
	if ( name_ == "C1c1" ) { get_symmops_C1c1(rt_out, cc ); }
	if ( name_ == "P12/m1" ) { get_symmops_P12slashm1(rt_out, cc ); }
	if ( name_ == "P121/m1" ) { get_symmops_P121slashm1(rt_out, cc ); }
	if ( name_ == "C12/m1" ) { get_symmops_C12slashm1(rt_out, cc ); }
	if ( name_ == "P12/c1" ) { get_symmops_P12slashc1(rt_out, cc ); }
	if ( name_ == "P121/c1" ) { get_symmops_P121slashc1(rt_out, cc ); }
	if ( name_ == "P121/n1" || name_ == "231" ) { get_symmops_P121slashn1(rt_out, cc ); }
	if ( name_ == "C12/c1" ) { get_symmops_C12slashc1(rt_out, cc ); }
	if ( name_ == "P222" ) { get_symmops_P222(rt_out, cc ); }
	if ( name_ == "P2221" ) { get_symmops_P2221(rt_out, cc ); }
	if ( name_ == "P21212" ) { get_symmops_P21212(rt_out, cc ); }
	if ( name_ == "P212121" ) { get_symmops_P212121(rt_out, cc ); }
	if ( name_ == "C2221" ) { get_symmops_C2221(rt_out, cc ); }
	if ( name_ == "C222" ) { get_symmops_C222(rt_out, cc ); }
	if ( name_ == "F222" ) { get_symmops_F222(rt_out, cc ); }
	if ( name_ == "I222" ) { get_symmops_I222(rt_out, cc ); }
	if ( name_ == "I212121" ) { get_symmops_I212121(rt_out, cc ); }
	if ( name_ == "Pmm2" ) { get_symmops_Pmm2(rt_out, cc ); }
	if ( name_ == "Pmc21" ) { get_symmops_Pmc21(rt_out, cc ); }
	if ( name_ == "Pcc2" ) { get_symmops_Pcc2(rt_out, cc ); }
	if ( name_ == "Pma2" ) { get_symmops_Pma2(rt_out, cc ); }
	if ( name_ == "Pca21" ) { get_symmops_Pca21(rt_out, cc ); }
	if ( name_ == "Pnc2" ) { get_symmops_Pnc2(rt_out, cc ); }
	if ( name_ == "Pmn21" ) { get_symmops_Pmn21(rt_out, cc ); }
	if ( name_ == "Pba2" ) { get_symmops_Pba2(rt_out, cc ); }
	if ( name_ == "Pna21" ) { get_symmops_Pna21(rt_out, cc ); }
	if ( name_ == "Pnn2" ) { get_symmops_Pnn2(rt_out, cc ); }
	if ( name_ == "Cmm2" ) { get_symmops_Cmm2(rt_out, cc ); }
	if ( name_ == "Cmc21" ) { get_symmops_Cmc21(rt_out, cc ); }
	if ( name_ == "Ccc2" ) { get_symmops_Ccc2(rt_out, cc ); }
	if ( name_ == "Amm2" ) { get_symmops_Amm2(rt_out, cc ); }
	if ( name_ == "Abm2" ) { get_symmops_Abm2(rt_out, cc ); }
	if ( name_ == "Ama2" ) { get_symmops_Ama2(rt_out, cc ); }
	if ( name_ == "Aba2" ) { get_symmops_Aba2(rt_out, cc ); }
	if ( name_ == "Fmm2" ) { get_symmops_Fmm2(rt_out, cc ); }
	if ( name_ == "Fdd2" ) { get_symmops_Fdd2(rt_out, cc ); }
	if ( name_ == "Imm2" ) { get_symmops_Imm2(rt_out, cc ); }
	if ( name_ == "Iba2" ) { get_symmops_Iba2(rt_out, cc ); }
	if ( name_ == "Ima2" ) { get_symmops_Ima2(rt_out, cc ); }
	if ( name_ == "Pmmm" ) { get_symmops_Pmmm(rt_out, cc ); }
	if ( name_ == "Pnnn:2" ) { get_symmops_Pnnn__2(rt_out, cc ); }
	if ( name_ == "Pccm" ) { get_symmops_Pccm(rt_out, cc ); }
	if ( name_ == "Pban:2" ) { get_symmops_Pban__2(rt_out, cc ); }
	if ( name_ == "Pmma" ) { get_symmops_Pmma(rt_out, cc ); }
	if ( name_ == "Pnna" ) { get_symmops_Pnna(rt_out, cc ); }
	if ( name_ == "Pmna" ) { get_symmops_Pmna(rt_out, cc ); }
	if ( name_ == "Pcca" ) { get_symmops_Pcca(rt_out, cc ); }
	if ( name_ == "Pbam" ) { get_symmops_Pbam(rt_out, cc ); }
	if ( name_ == "Pccn" ) { get_symmops_Pccn(rt_out, cc ); }
	if ( name_ == "Pbcm" ) { get_symmops_Pbcm(rt_out, cc ); }
	if ( name_ == "Pnnm" ) { get_symmops_Pnnm(rt_out, cc ); }
	if ( name_ == "Pmmn:2" ) { get_symmops_Pmmn__2(rt_out, cc ); }
	if ( name_ == "Pbcn" ) { get_symmops_Pbcn(rt_out, cc ); }
	if ( name_ == "Pbca" ) { get_symmops_Pbca(rt_out, cc ); }
	if ( name_ == "Pnma" ) { get_symmops_Pnma(rt_out, cc ); }
	if ( name_ == "Cmcm" ) { get_symmops_Cmcm(rt_out, cc ); }
	if ( name_ == "Cmca" ) { get_symmops_Cmca(rt_out, cc ); }
	if ( name_ == "Cmmm" ) { get_symmops_Cmmm(rt_out, cc ); }
	if ( name_ == "Cccm" ) { get_symmops_Cccm(rt_out, cc ); }
	if ( name_ == "Cmma" ) { get_symmops_Cmma(rt_out, cc ); }
	if ( name_ == "Ccca:2" ) { get_symmops_Ccca__2(rt_out, cc ); }
	if ( name_ == "Fmmm" ) { get_symmops_Fmmm(rt_out, cc ); }
	if ( name_ == "Fddd:2" ) { get_symmops_Fddd__2(rt_out, cc ); }
	if ( name_ == "Immm" ) { get_symmops_Immm(rt_out, cc ); }
	if ( name_ == "Ibam" ) { get_symmops_Ibam(rt_out, cc ); }
	if ( name_ == "Ibca" ) { get_symmops_Ibca(rt_out, cc ); }
	if ( name_ == "Imma" ) { get_symmops_Imma(rt_out, cc ); }
	if ( name_ == "P4" ) { get_symmops_P4(rt_out, cc ); }
	if ( name_ == "P41" ) { get_symmops_P41(rt_out, cc ); }
	if ( name_ == "P42" ) { get_symmops_P42(rt_out, cc ); }
	if ( name_ == "P43" ) { get_symmops_P43(rt_out, cc ); }
	if ( name_ == "I4" ) { get_symmops_I4(rt_out, cc ); }
	if ( name_ == "I41" ) { get_symmops_I41(rt_out, cc ); }
	if ( name_ == "P-4" ) { get_symmops_Pminus4(rt_out, cc ); }
	if ( name_ == "I-4" ) { get_symmops_Iminus4(rt_out, cc ); }
	if ( name_ == "P4/m" ) { get_symmops_P4slashm(rt_out, cc ); }
	if ( name_ == "P42/m" ) { get_symmops_P42slashm(rt_out, cc ); }
	if ( name_ == "P4/n:2" ) { get_symmops_P4slashn__2(rt_out, cc ); }
	if ( name_ == "P42/n:2" ) { get_symmops_P42slashn__2(rt_out, cc ); }
	if ( name_ == "I4/m" ) { get_symmops_I4slashm(rt_out, cc ); }
	if ( name_ == "I41/a:2" ) { get_symmops_I41slasha__2(rt_out, cc ); }
	if ( name_ == "P422" ) { get_symmops_P422(rt_out, cc ); }
	if ( name_ == "P4212" ) { get_symmops_P4212(rt_out, cc ); }
	if ( name_ == "P4122" ) { get_symmops_P4122(rt_out, cc ); }
	if ( name_ == "P41212" ) { get_symmops_P41212(rt_out, cc ); }
	if ( name_ == "P4222" ) { get_symmops_P4222(rt_out, cc ); }
	if ( name_ == "P42212" ) { get_symmops_P42212(rt_out, cc ); }
	if ( name_ == "P4322" ) { get_symmops_P4322(rt_out, cc ); }
	if ( name_ == "P43212" ) { get_symmops_P43212(rt_out, cc ); }
	if ( name_ == "I422" ) { get_symmops_I422(rt_out, cc ); }
	if ( name_ == "I4122" ) { get_symmops_I4122(rt_out, cc ); }
	if ( name_ == "P4mm" ) { get_symmops_P4mm(rt_out, cc ); }
	if ( name_ == "P4bm" ) { get_symmops_P4bm(rt_out, cc ); }
	if ( name_ == "P42cm" ) { get_symmops_P42cm(rt_out, cc ); }
	if ( name_ == "P42nm" ) { get_symmops_P42nm(rt_out, cc ); }
	if ( name_ == "P4cc" ) { get_symmops_P4cc(rt_out, cc ); }
	if ( name_ == "P4nc" ) { get_symmops_P4nc(rt_out, cc ); }
	if ( name_ == "P42mc" ) { get_symmops_P42mc(rt_out, cc ); }
	if ( name_ == "P42bc" ) { get_symmops_P42bc(rt_out, cc ); }
	if ( name_ == "I4mm" ) { get_symmops_I4mm(rt_out, cc ); }
	if ( name_ == "I4cm" ) { get_symmops_I4cm(rt_out, cc ); }
	if ( name_ == "I41md" ) { get_symmops_I41md(rt_out, cc ); }
	if ( name_ == "I41cd" ) { get_symmops_I41cd(rt_out, cc ); }
	if ( name_ == "P-42m" ) { get_symmops_Pminus42m(rt_out, cc ); }
	if ( name_ == "P-42c" ) { get_symmops_Pminus42c(rt_out, cc ); }
	if ( name_ == "P-421m" ) { get_symmops_Pminus421m(rt_out, cc ); }
	if ( name_ == "P-421c" ) { get_symmops_Pminus421c(rt_out, cc ); }
	if ( name_ == "P-4m2" ) { get_symmops_Pminus4m2(rt_out, cc ); }
	if ( name_ == "P-4c2" ) { get_symmops_Pminus4c2(rt_out, cc ); }
	if ( name_ == "P-4b2" ) { get_symmops_Pminus4b2(rt_out, cc ); }
	if ( name_ == "P-4n2" ) { get_symmops_Pminus4n2(rt_out, cc ); }
	if ( name_ == "I-4m2" ) { get_symmops_Iminus4m2(rt_out, cc ); }
	if ( name_ == "I-4c2" ) { get_symmops_Iminus4c2(rt_out, cc ); }
	if ( name_ == "I-42m" ) { get_symmops_Iminus42m(rt_out, cc ); }
	if ( name_ == "I-42d" ) { get_symmops_Iminus42d(rt_out, cc ); }
	if ( name_ == "P4/mmm" ) { get_symmops_P4slashmmm(rt_out, cc ); }
	if ( name_ == "P4/mcc" ) { get_symmops_P4slashmcc(rt_out, cc ); }
	if ( name_ == "P4/nbm:2" ) { get_symmops_P4slashnbm__2(rt_out, cc ); }
	if ( name_ == "P4/nnc:2" ) { get_symmops_P4slashnnc__2(rt_out, cc ); }
	if ( name_ == "P4/mbm" ) { get_symmops_P4slashmbm(rt_out, cc ); }
	if ( name_ == "P4/mnc" ) { get_symmops_P4slashmnc(rt_out, cc ); }
	if ( name_ == "P4/nmm:2" ) { get_symmops_P4slashnmm__2(rt_out, cc ); }
	if ( name_ == "P4/ncc:2" ) { get_symmops_P4slashncc__2(rt_out, cc ); }
	if ( name_ == "P42/mmc" ) { get_symmops_P42slashmmc(rt_out, cc ); }
	if ( name_ == "P42/mcm" ) { get_symmops_P42slashmcm(rt_out, cc ); }
	if ( name_ == "P42/nbc:2" ) { get_symmops_P42slashnbc__2(rt_out, cc ); }
	if ( name_ == "P42/nnm:2" ) { get_symmops_P42slashnnm__2(rt_out, cc ); }
	if ( name_ == "P42/mbc" ) { get_symmops_P42slashmbc(rt_out, cc ); }
	if ( name_ == "P42/mnm" ) { get_symmops_P42slashmnm(rt_out, cc ); }
	if ( name_ == "P42/nmc:2" ) { get_symmops_P42slashnmc__2(rt_out, cc ); }
	if ( name_ == "P42/ncm:2" ) { get_symmops_P42slashncm__2(rt_out, cc ); }
	if ( name_ == "I4/mmm" ) { get_symmops_I4slashmmm(rt_out, cc ); }
	if ( name_ == "I4/mcm" ) { get_symmops_I4slashmcm(rt_out, cc ); }
	if ( name_ == "I41/amd:2" ) { get_symmops_I41slashamd__2(rt_out, cc ); }
	if ( name_ == "I41/acd:2" ) { get_symmops_I41slashacd__2(rt_out, cc ); }
	if ( name_ == "P3" ) { get_symmops_P3(rt_out, cc ); }
	if ( name_ == "P31" ) { get_symmops_P31(rt_out, cc ); }
	if ( name_ == "P32" ) { get_symmops_P32(rt_out, cc ); }
	if ( name_ == "R3:H" || name_ == "H3" ) { get_symmops_R3__H(rt_out, cc ); }
	if ( name_ == "P-3" ) { get_symmops_Pminus3(rt_out, cc ); }
	if ( name_ == "R-3:H" ) { get_symmops_Rminus3__H(rt_out, cc ); }
	if ( name_ == "P312" ) { get_symmops_P312(rt_out, cc ); }
	if ( name_ == "P321" ) { get_symmops_P321(rt_out, cc ); }
	if ( name_ == "P3112" ) { get_symmops_P3112(rt_out, cc ); }
	if ( name_ == "P3121" ) { get_symmops_P3121(rt_out, cc ); }
	if ( name_ == "P3212" ) { get_symmops_P3212(rt_out, cc ); }
	if ( name_ == "P3221" ) { get_symmops_P3221(rt_out, cc ); }
	if ( name_ == "R32:H" || name_ == "H32" ) { get_symmops_R32__H(rt_out, cc ); }
	if ( name_ == "P3m1" ) { get_symmops_P3m1(rt_out, cc ); }
	if ( name_ == "P31m" ) { get_symmops_P31m(rt_out, cc ); }
	if ( name_ == "P3c1" ) { get_symmops_P3c1(rt_out, cc ); }
	if ( name_ == "P31c" ) { get_symmops_P31c(rt_out, cc ); }
	if ( name_ == "R3m:H" ) { get_symmops_R3m__H(rt_out, cc ); }
	if ( name_ == "R3c:H" ) { get_symmops_R3c__H(rt_out, cc ); }
	if ( name_ == "P-31m" ) { get_symmops_Pminus31m(rt_out, cc ); }
	if ( name_ == "P-31c" ) { get_symmops_Pminus31c(rt_out, cc ); }
	if ( name_ == "P-3m1" ) { get_symmops_Pminus3m1(rt_out, cc ); }
	if ( name_ == "P-3c1" ) { get_symmops_Pminus3c1(rt_out, cc ); }
	if ( name_ == "R-3m:H" ) { get_symmops_Rminus3m__H(rt_out, cc ); }
	if ( name_ == "R-3c:H" ) { get_symmops_Rminus3c__H(rt_out, cc ); }
	if ( name_ == "P6" ) { get_symmops_P6(rt_out, cc ); }
	if ( name_ == "P61" ) { get_symmops_P61(rt_out, cc ); }
	if ( name_ == "P65" ) { get_symmops_P65(rt_out, cc ); }
	if ( name_ == "P62" ) { get_symmops_P62(rt_out, cc ); }
	if ( name_ == "P64" ) { get_symmops_P64(rt_out, cc ); }
	if ( name_ == "P63" ) { get_symmops_P63(rt_out, cc ); }
	if ( name_ == "P-6" ) { get_symmops_Pminus6(rt_out, cc ); }
	if ( name_ == "P6/m" ) { get_symmops_P6slashm(rt_out, cc ); }
	if ( name_ == "P63/m" ) { get_symmops_P63slashm(rt_out, cc ); }
	if ( name_ == "P622" ) { get_symmops_P622(rt_out, cc ); }
	if ( name_ == "P6122" ) { get_symmops_P6122(rt_out, cc ); }
	if ( name_ == "P6522" ) { get_symmops_P6522(rt_out, cc ); }
	if ( name_ == "P6222" ) { get_symmops_P6222(rt_out, cc ); }
	if ( name_ == "P6422" ) { get_symmops_P6422(rt_out, cc ); }
	if ( name_ == "P6322" ) { get_symmops_P6322(rt_out, cc ); }
	if ( name_ == "P6mm" ) { get_symmops_P6mm(rt_out, cc ); }
	if ( name_ == "P6cc" ) { get_symmops_P6cc(rt_out, cc ); }
	if ( name_ == "P63cm" ) { get_symmops_P63cm(rt_out, cc ); }
	if ( name_ == "P63mc" ) { get_symmops_P63mc(rt_out, cc ); }
	if ( name_ == "P-6m2" ) { get_symmops_Pminus6m2(rt_out, cc ); }
	if ( name_ == "P-6c2" ) { get_symmops_Pminus6c2(rt_out, cc ); }
	if ( name_ == "P-62m" ) { get_symmops_Pminus62m(rt_out, cc ); }
	if ( name_ == "P-62c" ) { get_symmops_Pminus62c(rt_out, cc ); }
	if ( name_ == "P6/mmm" ) { get_symmops_P6slashmmm(rt_out, cc ); }
	if ( name_ == "P6/mcc" ) { get_symmops_P6slashmcc(rt_out, cc ); }
	if ( name_ == "P63/mcm" ) { get_symmops_P63slashmcm(rt_out, cc ); }
	if ( name_ == "P63/mmc" ) { get_symmops_P63slashmmc(rt_out, cc ); }
	if ( name_ == "P23" ) { get_symmops_P23(rt_out, cc ); }
	if ( name_ == "F23" ) { get_symmops_F23(rt_out, cc ); }
	if ( name_ == "I23" ) { get_symmops_I23(rt_out, cc ); }
	if ( name_ == "P213" ) { get_symmops_P213(rt_out, cc ); }
	if ( name_ == "I213" ) { get_symmops_I213(rt_out, cc ); }
	if ( name_ == "Pm-3" ) { get_symmops_Pmminus3(rt_out, cc ); }
	if ( name_ == "Pn-3:2" ) { get_symmops_Pnminus3__2(rt_out, cc ); }
	if ( name_ == "Fm-3" ) { get_symmops_Fmminus3(rt_out, cc ); }
	if ( name_ == "Fd-3:2" ) { get_symmops_Fdminus3__2(rt_out, cc ); }
	if ( name_ == "Im-3" ) { get_symmops_Imminus3(rt_out, cc ); }
	if ( name_ == "Pa-3" ) { get_symmops_Paminus3(rt_out, cc ); }
	if ( name_ == "Ia-3" ) { get_symmops_Iaminus3(rt_out, cc ); }
	if ( name_ == "P432" ) { get_symmops_P432(rt_out, cc ); }
	if ( name_ == "P4232" ) { get_symmops_P4232(rt_out, cc ); }
	if ( name_ == "F432" ) { get_symmops_F432(rt_out, cc ); }
	if ( name_ == "F4132" ) { get_symmops_F4132(rt_out, cc ); }
	if ( name_ == "I432" ) { get_symmops_I432(rt_out, cc ); }
	if ( name_ == "P4332" ) { get_symmops_P4332(rt_out, cc ); }
	if ( name_ == "P4132" ) { get_symmops_P4132(rt_out, cc ); }
	if ( name_ == "I4132" ) { get_symmops_I4132(rt_out, cc ); }
	if ( name_ == "P-43m" ) { get_symmops_Pminus43m(rt_out, cc ); }
	if ( name_ == "F-43m" ) { get_symmops_Fminus43m(rt_out, cc ); }
	if ( name_ == "I-43m" ) { get_symmops_Iminus43m(rt_out, cc ); }
	if ( name_ == "P-43n" ) { get_symmops_Pminus43n(rt_out, cc ); }
	if ( name_ == "F-43c" ) { get_symmops_Fminus43c(rt_out, cc ); }
	if ( name_ == "I-43d" ) { get_symmops_Iminus43d(rt_out, cc ); }
	if ( name_ == "Pm-3m" ) { get_symmops_Pmminus3m(rt_out, cc ); }
	if ( name_ == "Pn-3n:2" ) { get_symmops_Pnminus3n__2(rt_out, cc ); }
	if ( name_ == "Pm-3n" ) { get_symmops_Pmminus3n(rt_out, cc ); }
	if ( name_ == "Pn-3m:2" ) { get_symmops_Pnminus3m__2(rt_out, cc ); }
	if ( name_ == "Fm-3m" ) { get_symmops_Fmminus3m(rt_out, cc ); }
	if ( name_ == "Fm-3c" ) { get_symmops_Fmminus3c(rt_out, cc ); }
	if ( name_ == "Fd-3m:2" ) { get_symmops_Fdminus3m__2(rt_out, cc ); }
	if ( name_ == "Fd-3c:2" ) { get_symmops_Fdminus3c__2(rt_out, cc ); }
	if ( name_ == "Im-3m" ) { get_symmops_Imminus3m(rt_out, cc ); }
	if ( name_ == "Ia-3d" ) { get_symmops_Iaminus3d(rt_out, cc ); }
	if ( name_ == "B11m" ) { get_symmops_B11m(rt_out, cc ); }

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
	if ( name_ == "P121/n1" || name_ == "231" ) return "P 1 21/n 1";
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

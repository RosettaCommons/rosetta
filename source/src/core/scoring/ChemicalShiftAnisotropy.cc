// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ChemicalShiftAnisotropy.cc
/// @brief  Uses NMR CSA for scoring (Bertram R et al J of Magn Reson 147, 9-16)
/// @author Lei Shi

//Unit headers
#include <core/scoring/ChemicalShiftAnisotropy.hh>

// Package headers

// Project headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>

#include <basic/options/option.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

//C++ headers
#include <iostream>
#include <fstream>
#include <string>

/// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/csa.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

//Membrane code to find out membrane normal
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh>
#include <core/scoring/MembranePotential.hh>

static basic::Tracer tr( "core.scoring.ChemicalShiftAnisotropy" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

//////////////////////////////////////////////////////
//@brief reads in CSA data file
//////////////////////////////////////////////////////
extern void store_CSA_in_pose(ChemicalShiftAnisotropyOP csa_info, core::pose::Pose& pose) {
	pose.data().set(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA, csa_info);
}

extern ChemicalShiftAnisotropyCOP retrieve_CSA_from_pose( core::pose::Pose const& pose) {
	if ( pose.data().has(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) ) {
		return utility::pointer::static_pointer_cast< core::scoring::ChemicalShiftAnisotropy const > ( pose.data().get_const_ptr(
			core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) );
	};
	return nullptr;
}

extern ChemicalShiftAnisotropyOP retrieve_CSA_from_pose(core::pose::Pose& pose) {
	if ( pose.data().has(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) ) {
		return utility::pointer::static_pointer_cast< core::scoring::ChemicalShiftAnisotropy > ( pose.data().get_ptr(
			core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) );
	};
	return nullptr;
}

void ChemicalShiftAnisotropy::show(std::ostream& out) const {
	Size ct=0;
	for ( auto const & All_CSA_line : All_CSA_lines_ ) {
		out << "CSA "<<++ct << "     ";
		out << All_CSA_line << std::endl;
	}
}

void CSA::show(std::ostream& out) const {
	using namespace ObjexxFCL::format;
	out << RJ(4, res1_) << RJ(5, CSAval_);
}

std::ostream& operator<<(std::ostream& out, CSA const& csa) {
	csa.show(out);
	return out;
}

std::ostream& operator<<(std::ostream& out, ChemicalShiftAnisotropy const& csa) {
	csa.show(out);
	return out;
}

ChemicalShiftAnisotropy::ChemicalShiftAnisotropy(ChemicalShiftAnisotropy const& other) :
	basic::datacache::CacheableData(other) {
	All_CSA_lines_ = other.All_CSA_lines_;
}

//explicit assignment operator to initialize buffers
ChemicalShiftAnisotropy&
ChemicalShiftAnisotropy::operator=(ChemicalShiftAnisotropy const & other) {
	basic::datacache::CacheableData::operator=(other);
	All_CSA_lines_ = other.All_CSA_lines_;
	return *this;
}

void ChemicalShiftAnisotropy::read_CSA_file( std::string const& filename ) {

	std::string line;
	utility::io::izstream infile(filename.c_str());

	if ( !infile.good() ) {
		throw( utility::excn::EXCN_FileNotFound( filename ) );
	}

	//tr.Info << "Reading CSA file " << filename << std::endl;
	while ( getline(infile, line) ) {
		std::istringstream line_stream(line);
		std::string atom1;
		Size res1;
		Real sigma1, sigma2, sigma3, weight;
		Real CSAval,CSAerr;
		line_stream >> res1 >> atom1 >> sigma1 >> sigma2 >> sigma3 >> CSAval >> CSAerr;
		weight=3.0/(sigma1+sigma2+sigma3);

		if ( line_stream.fail() ) {
			tr.Error << "couldn't read line " << line << " in csa-file " << filename << std::endl;
			throw( utility::excn::EXCN_BadInput(" invalid line "+line+" in csa-file "+filename));
		}

		if ( atom1=="N" && sigma1<=sigma2 && sigma2<=sigma3 && res1>1 ) {
			All_CSA_lines_.push_back(CSA(res1,atom1,sigma1,sigma2,sigma3,CSAval,CSAerr,weight));
		} else {
			if ( atom1!="N" ) {
				throw( utility::excn::EXCN_BadInput( "only N15 CSA  is supported, not yet for " ));
			}
			if ( res1<=1 ) {
				throw( utility::excn::EXCN_BadInput( "residue id needs to be greater than 1" ));
			}
			if ( sigma1>sigma2 ||  sigma1>sigma2 ||  sigma2>sigma3 ) {
				throw( utility::excn::EXCN_BadInput( "sigma1 < sigma2 < sigma3" ));
			}
		}

	} //end of readline
}//end of read file

void ChemicalShiftAnisotropy::read_CSA_file( ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !option[ OptionKeys::in::file::csa ].user() ) {
		tr.Warning << "no CSA file specified" << std::endl;
		return;
	} else {
		std::string filename( option[ OptionKeys::in::file::csa ]()[1] );
		read_CSA_file( filename);
	}
}

//compute csascore
Real ChemicalShiftAnisotropy::compute_csascore(core::pose::Pose & pose) {
	Real total_dev=0.0;
	numeric::xyzVector<Real> memnorm;

	//Should we use the membrane normal?
	if ( basic::options::option[ basic::options::OptionKeys::csa::useZ].user() ) {
		memnorm = numeric::xyzVector<Real> (0.0,0.0,1.0);
		memnorm.normalize();
	} else {
		//add the membrane code here
		// core::scoring::MembraneEmbed membrane_embed( core::scoring::MembraneEmbed_from_pose( pose ));
		memnorm = core::scoring::nonconst_MembraneEmbed_from_pose( pose ).normal();
		memnorm.normalize();
	}
	tr.Trace << "memnorm.x(): " << memnorm.x() << " memnorm.y() " << memnorm.y() << " memnorm.z() " << memnorm.z() << std::endl;


	for ( auto & All_CSA_line : All_CSA_lines_ ) {
		//      tr.Trace << "it->res1(): " << it->res1() << " it->atom1() " << it->atom1() << std::endl;
		//      tr.Trace << "it->res2(): " << it->res2() << " it->atom2() " << it->atom2() << std::endl;
		//      tr.Trace << "it->res3(): " << it->res3() << " it->atom3() " << it->atom3() << std::endl;
		numeric::xyzVector<Real> v1( pose.residue(All_CSA_line.res2()).atom(All_CSA_line.atom2()).xyz() - pose.residue(All_CSA_line.res1()).atom(All_CSA_line.atom1()).xyz());
		numeric::xyzVector<Real> u1=v1.normalized();
		numeric::xyzVector<Real> v2( pose.residue(All_CSA_line.res3()).atom(All_CSA_line.atom3()).xyz() - pose.residue(All_CSA_line.res1()).atom(All_CSA_line.atom1()).xyz());
		numeric::xyzVector<Real> u2=v2.normalized();
		//the cross product order is flipped
		numeric::xyzVector<Real> b=(u2.cross_product(u1)).normalized();
		numeric::xyzVector<Real> n=b.cross_product(u1);
		numeric::xyzMatrix<Real> F(numeric::xyzMatrix<Real>::cols(u1.x(),u1.y(),u1.z(),n.x(),n.y(),n.z(),b.x(),b.y(),b.z()));

		//apply matrix rotation here
		numeric::xyzMatrix< Real > malpha = numeric::x_rotation_matrix_degrees( All_CSA_line.alpha() );
		numeric::xyzMatrix< Real > mbeta  = numeric::z_rotation_matrix_degrees( All_CSA_line.beta() );
		numeric::xyzMatrix< Real > PAF(F*mbeta*malpha);

		numeric::xyzMatrix<Real> sigmas(numeric::xyzMatrix<Real>::cols(All_CSA_line.sigma3(),0,0,0,All_CSA_line.sigma1(),0,0,0,All_CSA_line.sigma2()));


		//compute the observed CSA
		Real compCSA = All_CSA_line.CSAval_computed_=dot_product(memnorm,numeric::product(PAF*sigmas*PAF.transposed(),memnorm));

		//compute the difference and then energy
		Real devCSA = fabs(compCSA - All_CSA_line.CSAval()) < All_CSA_line.CSAerr() ? 0.0 : fabs(fabs(compCSA - All_CSA_line.CSAval())-All_CSA_line.CSAerr());
		//Add weight to energy and derivatives
		total_dev += All_CSA_line.weight()*devCSA*devCSA;

		tr.Trace << "resi: " << All_CSA_line.res1() << " compCSA: " << compCSA << " expCSA: " << All_CSA_line.CSAval() << " expErr " << All_CSA_line.CSAerr() << " devCSA: " << devCSA  << " total_dev "<< total_dev << std::endl;

		//derivatives
		numeric::xyzVector<Real> v0(numeric::xyzVector<Real> (0.0,0.0,0.0));
		numeric::xyzVector<Real> vx(numeric::xyzVector<Real> (1.0,0.0,0.0));
		numeric::xyzVector<Real> vy(numeric::xyzVector<Real> (0.0,1.0,0.0));
		numeric::xyzVector<Real> vz(numeric::xyzVector<Real> (0.0,0.0,1.0));

		//temporary data structure to compute derivatives
		numeric::xyzVector<Real> du1_dx(0.0);
		numeric::xyzVector<Real> du2_dx(0.0);
		Real p = 0.0;
		Real q = 0.0;
		Real r = 0.0;
		numeric::xyzMatrix<Real> S(0.0);
		numeric::xyzMatrix<Real> dF_dx(0.0);
		numeric::xyzMatrix<Real> dPAF_dx(0.0);
		Real dcompCSA_dx=0.0;

		//derivative values
		Real dE_xA=0.0;
		Real dE_yA=0.0;
		Real dE_zA=0.0;
		Real dE_xB=0.0;
		Real dE_yB=0.0;
		Real dE_zB=0.0;
		Real dE_xC=0.0;
		Real dE_yC=0.0;
		Real dE_zC=0.0;

		CSA& csa = All_CSA_line;
		if ( fabs(compCSA - All_CSA_line.CSAval()) <= All_CSA_line.CSAerr() ) {
			//    tr.Info << "deri 0: " << fabs(compCSA - it->CSAval()) << " , " << it->CSAerr() << std::endl;
			csa.f1ij_[0] = 0;
			csa.f1ij_[1] = 0;
			csa.f1ij_[2] = 0;
			csa.f2ij_[0] = 0;
			csa.f2ij_[1] = 0;
			csa.f2ij_[2] = 0;
			csa.f3ij_[0] = 0;
			csa.f3ij_[1] = 0;
			csa.f3ij_[2] = 0;
		} else if ( compCSA - All_CSA_line.CSAval() > All_CSA_line.CSAerr() ) {
			//    tr.Info << "deri 1: " << compCSA - it->CSAval() << " , " << it->CSAerr() << std::endl;

			//xA
			//swith order of u1 and u2 in the derivation
			du1_dx=-vx/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,-vx);
			du2_dx=-vx/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,-vx);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_xA=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f1ij_[0]=dE_xA;
			//tr.Info << "csa.f1ij_[0]: " << csa.f1ij_[0] << std::endl;

			//yA
			//swith order of u1 and u2 in the derivation
			du1_dx=-vy/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,-vy);
			du2_dx=-vy/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,-vy);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_yA=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f1ij_[1]=dE_yA;
			//tr.Info << "csa.f1ij_[1]: " << csa.f1ij_[1] << std::endl;

			//zA
			//swith order of u1 and u2 in the derivation
			du1_dx=-vz/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,-vz);
			du2_dx=-vz/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,-vz);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_zA=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f1ij_[2]=dE_zA;
			//tr.Info << "csa.f1ij_[2]: " << csa.f1ij_[2] << std::endl;

			//xB
			//swith order of u1 and u2 in the derivation
			du1_dx=vx/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,vx);
			du2_dx=v0/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,v0);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_xB=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f2ij_[0]=dE_xB;
			//tr.Info << "csa.f2ij_[0]: " << csa.f2ij_[0] << std::endl;

			//yB
			//swith order of u1 and u2 in the derivation
			du1_dx=vy/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,vy);
			du2_dx=v0/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,v0);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_yB=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f2ij_[1]=dE_yB;
			//tr.Info << "csa.f2ij_[1]: " << csa.f2ij_[1] << std::endl;

			//zB
			//swith order of u1 and u2 in the derivation
			du1_dx=vz/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,vz);
			du2_dx=v0/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,v0);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_zB=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f2ij_[2]=dE_zB;
			//tr.Info << "csa.f2ij_[2]: " << csa.f2ij_[2] << std::endl;

			//xC
			//swith order of u1 and u2 in the derivation
			du1_dx=v0/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,v0);
			du2_dx=vx/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,vx);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_xC=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f3ij_[0]=dE_xC;
			//tr.Info << "csa.f3ij_[0]: " << csa.f3ij_[0] << std::endl;

			//yC
			//swith order of u1 and u2 in the derivation
			du1_dx=v0/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,v0);
			du2_dx=vy/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,vy);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_yC=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f3ij_[1]=dE_yC;
			//tr.Info << "csa.f3ij_[1]: " << csa.f3ij_[1] << std::endl;

			//zC
			//swith order of u1 and u2 in the derivation
			du1_dx=v0/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,v0);
			du2_dx=vz/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,vz);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_zC=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() - All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f3ij_[2]=dE_zC;
			//tr.Info << "csa.f3ij_[2]: " << csa.f3ij_[2] << std::endl;

		} else {
			//tr.Info << "deri 2: " << compCSA - it->CSAval() << " , " << it->CSAerr() << std::endl;

			//xA
			//swith order of u1 and u2 in the derivation
			du1_dx=-vx/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,-vx);
			du2_dx=-vx/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,-vx);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_xA=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f1ij_[0]=dE_xA;
			//tr.Info << "csa.f1ij_[0]: " << csa.f1ij_[0] << std::endl;

			//yA
			//swith order of u1 and u2 in the derivation
			du1_dx=-vy/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,-vy);
			du2_dx=-vy/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,-vy);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_yA=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f1ij_[1]=dE_yA;
			//tr.Info << "csa.f1ij_[1]: " << csa.f1ij_[1] << std::endl;

			//zA
			//swith order of u1 and u2 in the derivation
			du1_dx=-vz/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,-vz);
			du2_dx=-vz/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,-vz);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_zA=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f1ij_[2]=dE_zA;
			//tr.Info << "csa.f1ij_[2]: " << csa.f1ij_[2] << std::endl;

			//xB
			//swith order of u1 and u2 in the derivation
			du1_dx=vx/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,vx);
			du2_dx=v0/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,v0);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_xB=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f2ij_[0]=dE_xB;
			//tr.Info << "csa.f2ij_[0]: " << csa.f2ij_[0] << std::endl;

			//yB
			//swith order of u1 and u2 in the derivation
			du1_dx=vy/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,vy);
			du2_dx=v0/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,v0);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_yB=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f2ij_[1]=dE_yB;
			//tr.Info << "csa.f2ij_[1]: " << csa.f2ij_[1] << std::endl;

			//zB
			//swith order of u1 and u2 in the derivation
			du1_dx=vz/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,vz);
			du2_dx=v0/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,v0);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_zB=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f2ij_[2]=dE_zB;
			//tr.Info << "csa.f2ij_[2]: " << csa.f2ij_[2] << std::endl;

			//xC
			//swith order of u1 and u2 in the derivation
			du1_dx=v0/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,v0);
			du2_dx=vx/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,vx);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_xC=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f3ij_[0]=dE_xC;
			//tr.Info << "csa.f3ij_[0]: " << csa.f3ij_[0] << std::endl;

			//yC
			//swith order of u1 and u2 in the derivation
			du1_dx=v0/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,v0);
			du2_dx=vy/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,vy);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_yC=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f3ij_[1]=dE_yC;
			//tr.Info << "csa.f3ij_[1]: " << csa.f3ij_[1] << std::endl;

			//zC
			//swith order of u1 and u2 in the derivation
			du1_dx=v0/v1.length()-v1*(1/(v1.length()*v1.length()*v1.length()))*dot_product(v1,v0);
			du2_dx=vz/v2.length()-v2*(1/(v2.length()*v2.length()*v2.length()))*dot_product(v2,vz);
			p=dot_product(u1,cross_product(b,du1_dx));
			q=(1/cross_product(u2,u1).length())*dot_product(u1,cross_product(u2,du1_dx));
			r=(1/cross_product(u2,u1).length())*dot_product(n,cross_product(du2_dx,u1)+cross_product(u2,du1_dx));
			S=numeric::xyzMatrix<Real>::cols(0,-p,-q,p,0,-r,q,r,0);
			dF_dx=F*S;
			dPAF_dx=dF_dx*mbeta*malpha;
			dcompCSA_dx=2*dot_product(memnorm,numeric::product(dPAF_dx*sigmas*PAF.transposed(),memnorm));
			dE_zC=All_CSA_line.weight()*2*( compCSA - All_CSA_line.CSAval() + All_CSA_line.CSAerr() )*dcompCSA_dx;
			csa.f3ij_[2]=dE_zC;
			//tr.Info << "csa.f3ij_[2]: " << csa.f3ij_[2] << std::endl;
		}
	}
	//tr.Info << "All_CSA_lines_.size(): " << All_CSA_lines_.size() << std::endl;

	//return energy
	return  total_dev;
}//end of compute_csascore

} //namespace Scoring
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::CSA::save( Archive & arc ) const {
	arc( CEREAL_NVP( CSAval_computed_ ) ); // Real
	arc( CEREAL_NVP( f1ij_ ) ); // core::Vector
	arc( CEREAL_NVP( f2ij_ ) ); // core::Vector
	arc( CEREAL_NVP( f3ij_ ) ); // core::Vector
	arc( CEREAL_NVP( res1_ ) ); // Size
	arc( CEREAL_NVP( res2_ ) ); // Size
	arc( CEREAL_NVP( res3_ ) ); // Size
	arc( CEREAL_NVP( atom1_ ) ); // std::string
	arc( CEREAL_NVP( atom2_ ) ); // std::string
	arc( CEREAL_NVP( atom3_ ) ); // std::string
	arc( CEREAL_NVP( sigma1_ ) ); // Real
	arc( CEREAL_NVP( sigma2_ ) ); // Real
	arc( CEREAL_NVP( sigma3_ ) ); // Real
	arc( CEREAL_NVP( alpha_ ) ); // Real
	arc( CEREAL_NVP( beta_ ) ); // Real
	arc( CEREAL_NVP( gamma_ ) ); // Real
	arc( CEREAL_NVP( CSAval_ ) ); // Real
	arc( CEREAL_NVP( CSAerr_ ) ); // Real
	arc( CEREAL_NVP( weight_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::CSA::load( Archive & arc ) {
	arc( CSAval_computed_ ); // Real
	arc( f1ij_ ); // core::Vector
	arc( f2ij_ ); // core::Vector
	arc( f3ij_ ); // core::Vector
	arc( res1_ ); // Size
	arc( res2_ ); // Size
	arc( res3_ ); // Size
	arc( atom1_ ); // std::string
	arc( atom2_ ); // std::string
	arc( atom3_ ); // std::string
	arc( sigma1_ ); // Real
	arc( sigma2_ ); // Real
	arc( sigma3_ ); // Real
	arc( alpha_ ); // Real
	arc( beta_ ); // Real
	arc( gamma_ ); // Real
	arc( CSAval_ ); // Real
	arc( CSAerr_ ); // Real
	arc( weight_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::CSA );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::ChemicalShiftAnisotropy::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( All_CSA_lines_ ) ); // CSA_lines
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::ChemicalShiftAnisotropy::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( All_CSA_lines_ ); // CSA_lines
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::ChemicalShiftAnisotropy );
CEREAL_REGISTER_TYPE( core::scoring::ChemicalShiftAnisotropy )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_ChemicalShiftAnisotropy )
#endif // SERIALIZATION



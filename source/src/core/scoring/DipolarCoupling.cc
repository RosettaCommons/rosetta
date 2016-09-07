// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/DipolarCoupling.cc
/// @brief  Uses NMR DC for scoring (Bertram R et al J of Magn Reson 147, 9-16)
/// @author Lei Shi

//Unit headers
#include <core/scoring/DipolarCoupling.hh>

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
#include <basic/options/keys/dc.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

//membrane
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh>
#include <core/scoring/MembranePotential.hh>

static THREAD_LOCAL basic::Tracer tr( "core.scoring.DipolarCoupling" );

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
//@brief reads in DC data file
//////////////////////////////////////////////////////
extern void store_DC_in_pose(DipolarCouplingOP dc_info, core::pose::Pose& pose) {
	pose.data().set(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA, dc_info);
}

extern DipolarCouplingCOP retrieve_DC_from_pose( core::pose::Pose const& pose) {
	if ( pose.data().has(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) ) {
		return utility::pointer::static_pointer_cast< core::scoring::DipolarCoupling const > ( pose.data().get_const_ptr(
			core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) );
	};
	return nullptr;
}

extern DipolarCouplingOP retrieve_DC_from_pose(core::pose::Pose& pose) {
	if ( pose.data().has(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) ) {
		return utility::pointer::static_pointer_cast< core::scoring::DipolarCoupling > ( pose.data().get_ptr(
			core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) );
	};
	return nullptr;
}

void DipolarCoupling::show(std::ostream& out) const {
	Size ct=0;
	for (const auto & All_DC_line : All_DC_lines_) {
		out << "DC " << ++ct << "     ";
		out << All_DC_line << std::endl;
	}
}

void DC::show(std::ostream& out) const {
	using namespace ObjexxFCL::format;
	out << RJ(4, res1_) << RJ(5, DCval_);
}

std::ostream& operator<<(std::ostream& out, DC const& dc) {
	dc.show(out);
	return out;
}

std::ostream& operator<<(std::ostream& out, DipolarCoupling const& dc) {
	dc.show(out);
	return out;
}

DipolarCoupling::DipolarCoupling(DipolarCoupling const& other) :
	basic::datacache::CacheableData(other) {
	All_DC_lines_ = other.All_DC_lines_;
}

//explicit asnumeric::signment operator to initialize buffers
DipolarCoupling&
DipolarCoupling::operator=(DipolarCoupling const & other) {
	basic::datacache::CacheableData::operator=(other);
	All_DC_lines_ = other.All_DC_lines_;
	return *this;
}

void DipolarCoupling::read_DC_file( std::string const& filename ) {

	std::string line;
	utility::io::izstream infile(filename.c_str());

	if ( !infile.good() ) {
		throw( utility::excn::EXCN_FileNotFound( filename ) );
	}

	tr.Info << "Reading DC file " << filename << std::endl;
	while ( getline(infile, line) ) {
		std::istringstream line_stream(line);
		std::string atom1,atom2;
		Size res1,res2;
		Real weight;
		Real DCval,DCerr;
		line_stream >> res1 >> atom1 >> res2 >> atom2 >> DCval >> DCerr;
		weight=1.0;

		if ( line_stream.fail() ) {
			tr.Error << "couldn't read line " << line << " in dc-file " << filename << std::endl;
			throw( utility::excn::EXCN_BadInput(" invalid line "+line+" in dc-file "+filename));
		} else {
			All_DC_lines_.push_back(DC(res1,atom1,res2,atom2,DCval,DCerr,weight));
		}
	} //end of readline
}//end of read file

void DipolarCoupling::read_DC_file( ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !option[ OptionKeys::in::file::dc ].user() ) {
		tr.Warning << "no DC file specified" << std::endl;
		return;
	} else {
		std::string filename( option[ OptionKeys::in::file::dc ]()[1] );
		read_DC_file( filename);
	}
}


std::string element_string_dc(std::string atom) {
	if ( atom == "HN" || atom == "H" || atom == "HA" ) {
		return "H";
	}
	if ( atom == "C" || atom == "CA" ) {
		return "C";
	}
	if ( atom == "N" ) {
		return "N";
	}
	throw(utility::excn::EXCN_BadInput("unknown atom for DC: " + atom));
	return ""; //to make compile happy.
}
//////////////////////////////////////////////////////
//@brief returns type of DC data N-H, HA-CA etc
//////////////////////////////////////////////////////
DC::DC_TYPE DC::get_DC_data_type(std::string const & atom1, std::string const & atom2) {
	std::string elem1(element_string_dc(atom1));
	std::string elem2(element_string_dc(atom2));

	DC_TYPE DC_type;
	if ( (elem1 == "N" && elem2 == "H") || (elem1 == "H" && elem2 == "N") ) {
		DC_type = DC_TYPE_NH;
	} else if ( (elem1 == "C" && elem2 == "H") || (elem1 == "H" && elem2 == "C") ) {
		DC_type = DC_TYPE_CH;
	} else if ( (elem1 == "C" && elem2 == "N") || (elem1 == "N" && elem2 == "C") ) {
		DC_type = DC_TYPE_NC;
	} else if ( (elem1 == "C" && elem2 == "C") ) {
		DC_type = DC_TYPE_CC;
	} else {
		throw(utility::excn::EXCN_BadInput(
			"unknown combination of atoms for DC " + atom1 + " " + atom2));
	}
	return DC_type;
}

//compute dcscore
Real DipolarCoupling::compute_dcscore(core::pose::Pose & pose) {
	Real total_dev=0.0;
	numeric::xyzVector<Real> memnorm;
	//Should we use the membrane normal?
	if ( basic::options::option[ basic::options::OptionKeys::dc::useZ].user() ) {
		memnorm=numeric::xyzVector<Real> (0.0,0.0,1.0);
		memnorm.normalize();
	} else {
		//add the membrane code here!!
		memnorm = core::scoring::nonconst_MembraneEmbed_from_pose( pose ).normal();
		memnorm.normalize();
	}
	tr.Trace << "memnorm.x(): " << memnorm.x() << " memnorm.y() " << memnorm.y() << " memnorm.z() " << memnorm.z() << std::endl;

	for (auto & All_DC_line : All_DC_lines_) {
		//        tr.Trace << "it->res1(): " << it->res1() << " it->atom1() " << it->atom1() << std::endl;
		//        tr.Trace << "it->res2(): " << it->res2() << " it->atom2() " << it->atom2() << std::endl;
		numeric::xyzVector<Real> v1( pose.residue(All_DC_line.res2()).atom(All_DC_line.atom2()).xyz() - pose.residue(All_DC_line.res1()).atom(All_DC_line.atom1()).xyz());
		numeric::xyzVector<Real> u1=v1.normalized();

		Real costheta=dot_product(u1,memnorm);

		//compute the observed DC
		Real compDC=All_DC_line.DCval_computed_=All_DC_line.Dconst()*(3*costheta*costheta-1);

		//compute the difference and then energy
		Real devDC;
		Real dcos_dxA=0.0;
		Real dcos_dyA=0.0;
		Real dcos_dzA=0.0;
		Real dcos_dxB=0.0;
		Real dcos_dyB=0.0;
		Real dcos_dzB=0.0;

		Real dv_xA=0.0;
		Real dv_yA=0.0;
		Real dv_zA=0.0;
		Real dv_xB=0.0;
		Real dv_yB=0.0;
		Real dv_zB=0.0;

		DC& dc = All_DC_line;

		if ( All_DC_line.DCval()<=All_DC_line.Dconst() ) {
			//devDC = fabs(fabs(compDC) - it->DCval()) < it->DCerr() ? 0.0 : fabs(fabs(compDC - it->DCval())-it->DCerr());
			//derivatives
			if ( fabs(fabs(compDC) - All_DC_line.DCval()) <= All_DC_line.DCerr() ) {
				devDC=0.0;
				dc.f1ij_[0] = 0;
				dc.f1ij_[1] = 0;
				dc.f1ij_[2] = 0;
				dc.f2ij_[0] = 0;
				dc.f2ij_[1] = 0;
				dc.f2ij_[2] = 0;
			} else {

				if ( fabs(compDC) - All_DC_line.DCval()>All_DC_line.DCerr() ) {
					devDC=fabs(compDC) - All_DC_line.DCval() - All_DC_line.DCerr();
				} else {
					devDC=fabs(compDC) - All_DC_line.DCval() + All_DC_line.DCerr();
				}
				dcos_dxA= ( -memnorm.x()/v1.length()+dot_product(v1,memnorm)*v1.x()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dyA= ( -memnorm.y()/v1.length()+dot_product(v1,memnorm)*v1.y()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dzA= ( -memnorm.z()/v1.length()+dot_product(v1,memnorm)*v1.z()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dxB= ( memnorm.x()/v1.length()-dot_product(v1,memnorm)*v1.x()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dyB= ( memnorm.y()/v1.length()-dot_product(v1,memnorm)*v1.y()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dzB= ( memnorm.z()/v1.length()-dot_product(v1,memnorm)*v1.z()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();

				dv_xA=6*All_DC_line.Dconst()*costheta*dcos_dxA;
				dv_yA=6*All_DC_line.Dconst()*costheta*dcos_dyA;
				dv_zA=6*All_DC_line.Dconst()*costheta*dcos_dzA;
				dv_xB=6*All_DC_line.Dconst()*costheta*dcos_dxB;
				dv_yB=6*All_DC_line.Dconst()*costheta*dcos_dyB;
				dv_zB=6*All_DC_line.Dconst()*costheta*dcos_dzB;

				dc.f1ij_[0] = All_DC_line.weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_xA;
				dc.f1ij_[1] = All_DC_line.weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_yA;
				dc.f1ij_[2] = All_DC_line.weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_zA;
				dc.f2ij_[0] = All_DC_line.weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_xB;
				dc.f2ij_[1] = All_DC_line.weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_yB;
				dc.f2ij_[2] = All_DC_line.weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_zB;
			}
		} else if ( All_DC_line.DCval()>All_DC_line.Dconst() && All_DC_line.DCval()<=2*All_DC_line.Dconst() ) {
			if ( fabs(compDC- All_DC_line.DCval()) <= All_DC_line.DCerr() ) {
				devDC=0.0;
				dc.f1ij_[0] = 0;
				dc.f1ij_[1] = 0;
				dc.f1ij_[2] = 0;
				dc.f2ij_[0] = 0;
				dc.f2ij_[1] = 0;
				dc.f2ij_[2] = 0;
			} else {

				if ( compDC- All_DC_line.DCval()>All_DC_line.DCerr() ) {
					devDC=compDC- All_DC_line.DCval() - All_DC_line.DCerr();
				} else {
					devDC=compDC- All_DC_line.DCval() + All_DC_line.DCerr();
				}

				dcos_dxA= ( -memnorm.x()/v1.length()+dot_product(v1,memnorm)*v1.x()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dyA= ( -memnorm.y()/v1.length()+dot_product(v1,memnorm)*v1.y()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dzA= ( -memnorm.z()/v1.length()+dot_product(v1,memnorm)*v1.z()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dxB= ( memnorm.x()/v1.length()-dot_product(v1,memnorm)*v1.x()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dyB= ( memnorm.y()/v1.length()-dot_product(v1,memnorm)*v1.y()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dzB= ( memnorm.z()/v1.length()-dot_product(v1,memnorm)*v1.z()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();

				dv_xA=6*All_DC_line.Dconst()*costheta*dcos_dxA;
				dv_yA=6*All_DC_line.Dconst()*costheta*dcos_dyA;
				dv_zA=6*All_DC_line.Dconst()*costheta*dcos_dzA;
				dv_xB=6*All_DC_line.Dconst()*costheta*dcos_dxB;
				dv_yB=6*All_DC_line.Dconst()*costheta*dcos_dyB;
				dv_zB=6*All_DC_line.Dconst()*costheta*dcos_dzB;

				dc.f1ij_[0] = All_DC_line.weight()*2*devDC*dv_xA;
				dc.f1ij_[1] = All_DC_line.weight()*2*devDC*dv_yA;
				dc.f1ij_[2] = All_DC_line.weight()*2*devDC*dv_zA;
				dc.f2ij_[0] = All_DC_line.weight()*2*devDC*dv_xB;
				dc.f2ij_[1] = All_DC_line.weight()*2*devDC*dv_yB;
				dc.f2ij_[2] = All_DC_line.weight()*2*devDC*dv_zB;
			}
		} else {
			throw(utility::excn::EXCN_BadInput( "DC value not in resonable range"));
		}

		//Add weight to energy
		total_dev+=All_DC_line.weight()*devDC*devDC;

		if ( tr.Trace.visible() ) {
			tr.Trace << "resi: " << All_DC_line.res1() << " compDC: " << compDC << " expDC: " << All_DC_line.DCval() << " expErr " << All_DC_line.DCerr() << " devDC: " << devDC << " total_dev "<< total_dev << std::endl;
		}
	}

	//return energy
	return  total_dev;
}//end of compute_dcscore

} //namespace Scoring
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::DC::save( Archive & arc ) const {
	arc( CEREAL_NVP( DCval_computed_ ) ); // Real
	arc( CEREAL_NVP( f1ij_ ) ); // core::Vector
	arc( CEREAL_NVP( f2ij_ ) ); // core::Vector
	arc( CEREAL_NVP( type_ ) ); // enum core::scoring::DC::DC_TYPE
	arc( CEREAL_NVP( res1_ ) ); // Size
	arc( CEREAL_NVP( res2_ ) ); // Size
	arc( CEREAL_NVP( atom1_ ) ); // std::string
	arc( CEREAL_NVP( atom2_ ) ); // std::string
	arc( CEREAL_NVP( DCval_ ) ); // Real
	arc( CEREAL_NVP( DCerr_ ) ); // Real
	arc( CEREAL_NVP( weight_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::DC::load( Archive & arc ) {
	arc( DCval_computed_ ); // Real
	arc( f1ij_ ); // core::Vector
	arc( f2ij_ ); // core::Vector
	arc( type_ ); // enum core::scoring::DC::DC_TYPE
	arc( res1_ ); // Size
	arc( res2_ ); // Size
	arc( atom1_ ); // std::string
	arc( atom2_ ); // std::string
	arc( DCval_ ); // Real
	arc( DCerr_ ); // Real
	arc( weight_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::DC );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::DipolarCoupling::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( All_DC_lines_ ) ); // DC_lines
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::DipolarCoupling::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( All_DC_lines_ ); // DC_lines
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::DipolarCoupling );
CEREAL_REGISTER_TYPE( core::scoring::DipolarCoupling )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_DipolarCoupling )
#endif // SERIALIZATION



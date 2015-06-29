// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

static thread_local basic::Tracer tr( "core.scoring.DipolarCoupling" );

namespace core {
namespace scoring {

//////////////////////////////////////////////////////
//@brief reads in DC data file
//////////////////////////////////////////////////////
extern void store_DC_in_pose(DipolarCouplingOP dc_info, core::pose::Pose& pose) {
	pose.data().set(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA, dc_info);
}

extern DipolarCouplingCOP retrieve_DC_from_pose( core::pose::Pose const& pose) {
	if (pose.data().has(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA)) {
		return utility::pointer::static_pointer_cast< core::scoring::DipolarCoupling const > ( pose.data().get_const_ptr(
				core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) );
	};
	return NULL;
}

extern DipolarCouplingOP retrieve_DC_from_pose(core::pose::Pose& pose) {
	if (pose.data().has(core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA)) {
		return utility::pointer::static_pointer_cast< core::scoring::DipolarCoupling > ( pose.data().get_ptr(
				core::pose::datacache::CacheableDataType::CHEMICAL_SHIFT_ANISOTROPY_DATA) );
	};
	return NULL;
}

void DipolarCoupling::show(std::ostream& out) const {
  Size ct=0;
  for (DC_lines::const_iterator it=All_DC_lines_.begin();it!=All_DC_lines_.end();it++){
        out << "DC "<<++ct << "     ";
    out << (*it) << std::endl;
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
	while (getline(infile, line)) {
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
        if (atom == "HN" || atom == "H" || atom == "HA")
                return "H";
        if (atom == "C" || atom == "CA")
                return "C";
        if (atom == "N")
                return "N";
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
        if ((elem1 == "N" && elem2 == "H") || (elem1 == "H" && elem2 == "N"))
                DC_type = DC_TYPE_NH;
        else if ((elem1 == "C" && elem2 == "H") || (elem1 == "H" && elem2 == "C"))
                DC_type = DC_TYPE_CH;
        else if ((elem1 == "C" && elem2 == "N") || (elem1 == "N" && elem2 == "C"))
                DC_type = DC_TYPE_NC;
        else if ((elem1 == "C" && elem2 == "C"))
                DC_type = DC_TYPE_CC;
        else
                throw(utility::excn::EXCN_BadInput(
                                "unknown combination of atoms for DC " + atom1 + " " + atom2));
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

	for (utility::vector1<core::scoring::DC>::iterator it = All_DC_lines_.begin(); it != All_DC_lines_.end(); ++it) {
//        tr.Trace << "it->res1(): " << it->res1() << " it->atom1() " << it->atom1() << std::endl;
//        tr.Trace << "it->res2(): " << it->res2() << " it->atom2() " << it->atom2() << std::endl;
			numeric::xyzVector<Real> v1( pose.residue(it->res2()).atom(it->atom2()).xyz() - pose.residue(it->res1()).atom(it->atom1()).xyz());
			numeric::xyzVector<Real> u1=v1.normalized();
			
			Real costheta=dot_product(u1,memnorm);

		//compute the observed DC
			Real compDC=it->DCval_computed_=it->Dconst()*(3*costheta*costheta-1);

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

      DC& dc = *it;

			if ( it->DCval()<=it->Dconst() ) {
      				//devDC = fabs(fabs(compDC) - it->DCval()) < it->DCerr() ? 0.0 : fabs(fabs(compDC - it->DCval())-it->DCerr());
    				//derivatives
				if (fabs(fabs(compDC) - it->DCval()) <= it->DCerr()) {
					devDC=0.0;
					dc.f1ij_[0] = 0;
					dc.f1ij_[1] = 0;
					dc.f1ij_[2] = 0;
					dc.f2ij_[0] = 0;
					dc.f2ij_[1] = 0;
					dc.f2ij_[2] = 0;
				} else {
					
					if ( fabs(compDC) - it->DCval()>it->DCerr() ) {
						devDC=fabs(compDC) - it->DCval() - it->DCerr();
					} else {
						devDC=fabs(compDC) - it->DCval() + it->DCerr();
					}
					dcos_dxA= ( -memnorm.x()/v1.length()+dot_product(v1,memnorm)*v1.x()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
					dcos_dyA= ( -memnorm.y()/v1.length()+dot_product(v1,memnorm)*v1.y()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
					dcos_dzA= ( -memnorm.z()/v1.length()+dot_product(v1,memnorm)*v1.z()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
					dcos_dxB= ( memnorm.x()/v1.length()-dot_product(v1,memnorm)*v1.x()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
					dcos_dyB= ( memnorm.y()/v1.length()-dot_product(v1,memnorm)*v1.y()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
					dcos_dzB= ( memnorm.z()/v1.length()-dot_product(v1,memnorm)*v1.z()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();

					dv_xA=6*it->Dconst()*costheta*dcos_dxA;
					dv_yA=6*it->Dconst()*costheta*dcos_dyA;
					dv_zA=6*it->Dconst()*costheta*dcos_dzA;
					dv_xB=6*it->Dconst()*costheta*dcos_dxB;
					dv_yB=6*it->Dconst()*costheta*dcos_dyB;
					dv_zB=6*it->Dconst()*costheta*dcos_dzB;

					dc.f1ij_[0] = it->weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_xA;
					dc.f1ij_[1] = it->weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_yA;
					dc.f1ij_[2] = it->weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_zA;
					dc.f2ij_[0] = it->weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_xB;
					dc.f2ij_[1] = it->weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_yB;
					dc.f2ij_[2] = it->weight()*numeric::sign(3*costheta*costheta-1)*2*devDC*dv_zB;
				}
			}
			else if ( it->DCval()>it->Dconst() && it->DCval()<=2*it->Dconst() ) {
				if (fabs(compDC- it->DCval()) <= it->DCerr()) {
					devDC=0.0;
					dc.f1ij_[0] = 0;
					dc.f1ij_[1] = 0;
					dc.f1ij_[2] = 0;
					dc.f2ij_[0] = 0;
					dc.f2ij_[1] = 0;
					dc.f2ij_[2] = 0;
				} else {

				if ( compDC- it->DCval()>it->DCerr() ) {
					devDC=compDC- it->DCval() - it->DCerr();
				} else {
					devDC=compDC- it->DCval() + it->DCerr();
				}

				dcos_dxA= ( -memnorm.x()/v1.length()+dot_product(v1,memnorm)*v1.x()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dyA= ( -memnorm.y()/v1.length()+dot_product(v1,memnorm)*v1.y()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dzA= ( -memnorm.z()/v1.length()+dot_product(v1,memnorm)*v1.z()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dxB= ( memnorm.x()/v1.length()-dot_product(v1,memnorm)*v1.x()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dyB= ( memnorm.y()/v1.length()-dot_product(v1,memnorm)*v1.y()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();
				dcos_dzB= ( memnorm.z()/v1.length()-dot_product(v1,memnorm)*v1.z()/(v1.length()*v1.length()*v1.length()) ) / memnorm.length();

				dv_xA=6*it->Dconst()*costheta*dcos_dxA;
				dv_yA=6*it->Dconst()*costheta*dcos_dyA;
				dv_zA=6*it->Dconst()*costheta*dcos_dzA;
				dv_xB=6*it->Dconst()*costheta*dcos_dxB;
				dv_yB=6*it->Dconst()*costheta*dcos_dyB;
				dv_zB=6*it->Dconst()*costheta*dcos_dzB;

				dc.f1ij_[0] = it->weight()*2*devDC*dv_xA;
				dc.f1ij_[1] = it->weight()*2*devDC*dv_yA;
				dc.f1ij_[2] = it->weight()*2*devDC*dv_zA;
				dc.f2ij_[0] = it->weight()*2*devDC*dv_xB;
				dc.f2ij_[1] = it->weight()*2*devDC*dv_yB;
				dc.f2ij_[2] = it->weight()*2*devDC*dv_zB;
			}
		} else
			throw(utility::excn::EXCN_BadInput( "DC value not in resonable range"));

		//Add weight to energy
			total_dev+=it->weight()*devDC*devDC;

			if ( tr.Trace.visible() ) {
			tr.Trace << "resi: " << it->res1() << " compDC: " << compDC << " expDC: " << it->DCval() << " expErr " << it->DCerr() << " devDC: " << devDC << " total_dev "<< total_dev << std::endl;
		}
	}

	//return energy
	return  total_dev;
	//return  total_dev/All_DC_lines_.size();

}//end of compute_dcscore

} //namespace Scoring
} //namespace core

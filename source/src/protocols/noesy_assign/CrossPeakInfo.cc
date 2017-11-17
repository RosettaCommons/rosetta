// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/CrossPeakInfo.hh>

// Package Headers
#include <protocols/noesy_assign/Exceptions.hh>

// Project Headers
//#include <core/chemical/AA.hh>

// Utility headers

// #include <utility/excn/Exceptions.hh>
// #include <utility/vector1.fwd.hh>
// #include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>
// Third-party Headers
#include <boost/functional/hash.hpp>


static basic::Tracer tr( "protocols.noesy_assign.crosspeaks" );
static basic::Tracer tr_labels( "protocols.noesy_assign.crosspeaks.labels" );
using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {

/// @details Auto-generated virtual destructor
CrossPeakInfo::~CrossPeakInfo() = default;

void CrossPeakInfo::show( std::ostream& os ) const {
	os << "CROSSPEAK: "
		<< proton_atom_name_ << " " << label_atom_type_ << " "
		<< "TOL: " << proton_tolerance_ << " " << label_tolerance_
		<< " from file " << filename_;
	os << "PROTON: ";
	fold_proton_resonance_.show( os );
	os << "LABEL: ";
	fold_label_resonance_.show( os );
	os << "MAX_NOE_DIST: " << max_noe_distance_;
}

void CrossPeakInfo::set_filename( std::string filename ) {
	filename_ = filename;
	exp_hash_ = boost::hash_value( filename );
}


std::ostream& operator<< ( std::ostream& os, CrossPeakInfo const& cpi ) {
	cpi.show( os );
	return os;
}


std::string CrossPeakInfo::label_atom_name( std::string const& proton_name, core::chemical::AA aa ) const {
	using namespace core::chemical; //for AA
	if ( label_atom_type_ == "N" || label_atom_type_ == "NC" ) {
		if ( aa == aa_arg ) {
			if ( proton_name == "HE" ) return "NE";
			if ( proton_name.substr(0,2) == "HH" ) return "N"+proton_name.substr(1,2); //HH11 HH12 HH21 HH22
		}
		if ( aa == aa_lys ) {
			if ( proton_name.substr(0,2) == "HZ" ) return "NZ";
		}
		if ( aa == aa_gln || aa == aa_his ) {
			if ( proton_name.substr(0,3) == "HE2" ) return "NE2"; //HE21, HE22
		}
		if ( aa == aa_asn ) {
			if ( proton_name.substr(0,3) == "HD2" ) return "ND2"; //HD21, HD22
		}
		if ( aa == aa_trp ) {
			if ( proton_name == "HE1" ) return "NE1";
		}
		if ( aa == aa_his ) {
			if ( proton_name == "HD1" ) return "ND1";
		}
		if ( proton_name == "H" ) return "N";
	} // atom type is "N"

	std::string name;
	if ( label_atom_type_ == "C" || label_atom_type_ == "NC" ) {
		if ( proton_name[ 0 ] == 'Q' && proton_name[ 1 ] == 'Q' ) { //QQX
			name = "C" + proton_name.substr(2,1);
			return name;
		}
		if ( proton_name[ 0 ] == 'Q' ) { //Qxx
			name = proton_name;
			name[ 0 ] ='C';
			return name;
		} if ( proton_name.substr(1,2) == "HB" ) {
			return "CB";
		} if ( proton_name.substr(1,2) == "HA" ) {
			return "CA";
		}
		if ( aa == aa_trp ) {
			if ( proton_name == "HH2" ) return "CH2";
			if ( proton_name == "HZ2" ) return "CZ2";
			if ( proton_name == "HZ3" ) return "CZ3";
			if ( proton_name == "HE3" ) return "CE3";
			if ( proton_name == "HD1" ) return "CD1";
		}
		if ( aa == aa_phe || aa == aa_tyr || aa == aa_his ) {
			if ( proton_name == "HZ" ) return "CZ";
			if ( proton_name.substr(0,2) == "HD" || proton_name.substr(0,2) == "HE"  ) return "C"+proton_name.substr(1,2);
		}

		if ( proton_name.substr(0,2) == "HA" ) return "CA";
		if ( proton_name.substr(0,2) == "HB" ) return "CB";
		if ( aa != aa_asn ) { //don't make HD21 -> CD substition for ASN
			Size len=proton_name.size()-2;
			if ( proton_name.substr(0,2) == "HG" || proton_name.substr(0,2)=="HD" || proton_name.substr(0,2)=="HE" ) return "C"+proton_name.substr(1,len < 1 ? 1 : len );
		}
		// already as aa!=asn  if ( proton_name.substr(0,2)== "HE" ) return "CE";
		// already as aa!=asn  if ( proton_name.substr(0,2)== "HD" ) return "CD";
	}
	if ( tr_labels.Trace.visible() ) {
		tr_labels.Trace << "proton_name " + proton_name + " not recognized for " + label_atom_type_ + " label on " + name_from_aa( aa ) << std::endl;
	}
	throw CREATE_EXCEPTION(EXCN_UnknownAtomname, "");
	return "no_atom";
}

}
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
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
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/string_util.hh>
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
#include <iostream>
#include <cstdlib>
#include <string>


static basic::Tracer tr("protocols.noesy_assign.crosspeaks");

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {

std::string CrossPeakInfo::label_atom_name( std::string const& proton_name, core::chemical::AA aa ) const {
  using namespace core::chemical; //for AA
  if ( label_atom_type_ == "N" ) {
    if ( aa == aa_arg ) {
      if ( proton_name == "HE" ) return "NE";
      if ( proton_name.substr(0,2) == "HH" ) return "N"+proton_name.substr(1,2); //HH11 HH12 HH21 HH22
    }
    if ( aa == aa_lys ) {
      if ( proton_name.substr(0,2) == "HZ" ) return "NZ";
    }
    if ( aa == aa_gln ) {
      if ( proton_name.substr(0,3) == "HE2" ) return "NE2"; //HE21, HE22
    }
    if ( aa == aa_asn ) {
      if ( proton_name.substr(0,3) == "HD2" ) return "ND2"; //HD21, HD22
    }
    if ( aa == aa_trp ) {
      if ( proton_name == "HE1" ) return "NE1";
    }
    if ( proton_name == "H" ) return "N";
  } // atom type is "N"

  std::string name;
  if ( label_atom_type_ == "C" ) {
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
    if ( aa == aa_phe || aa == aa_tyr ) {
      if ( proton_name == "HZ" ) return "CZ";
      if ( proton_name.substr(0,2) == "HD" || proton_name.substr(0,2) == "HE"  ) return "C"+proton_name.substr(1,2);
    }

    if ( proton_name.substr(0,2) == "HA" ) return "CA";
    if ( proton_name.substr(0,2) == "HB" ) return "CB";
    if ( aa != aa_asn ) { //don't make HD21 -> CD substition for ASN
			Size len=proton_name.size()-2;
      if ( proton_name.substr(0,2) == "HG" || proton_name.substr(0,2)=="HD" ) return "C"+proton_name.substr(1,len < 1 ? 1 : len );
    }


  }
  throw EXCN_UnknownAtomname("proton_name " + proton_name + " not recognized for " + label_atom_type_ + " label" );
  return "no_atom";
}

}
}

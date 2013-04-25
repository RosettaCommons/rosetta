// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NativeResidue.cc
///
/// @brief

/// @author Tim jacobs

//Core
#include <core/types.hh>

//Unit Headers
#include <devel/helixAssembly/NativeResidue.hh>

//Utility
#include <utility/string_util.hh>

//Devel
#include <devel/helixAssembly/NativeAtom.hh>

NativeResidue::NativeResidue():
res_type_(),
atoms_()
{}

NativeResidue::NativeResidue(const std::string & res_type, const std::vector<NativeAtom> & atoms):
res_type_(res_type),
atoms_(atoms)
{}

NativeResidue::~NativeResidue(){}

std::string NativeResidue::res_type(){
  return res_type_;
}

std::string NativeResidue::print() const{
  std::string output = "RESIDUE " + res_type_ + "\n";

  for(core::Size i=0; i<atoms_.size(); ++i){
      output += atoms_[i].print();
  }
  return output;
}

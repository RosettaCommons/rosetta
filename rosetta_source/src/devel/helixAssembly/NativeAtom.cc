// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NativeAtom.cc
///
/// @brief

/// @author Tim Jacobs

//Unit headers
#include <devel/helixAssembly/NativeAtom.hh>

//Utility
#include <utility/string_util.hh>

NativeAtom::NativeAtom():
type_(0),
x_(0),
y_(0),
z_(0)
{}

NativeAtom::NativeAtom(const unsigned short type, const float x, const float y, const float z):
type_(type),
x_(x),
y_(y),
z_(z)
{}

NativeAtom::~NativeAtom(){}

float NativeAtom::x(){
  return x_;
}

float NativeAtom::y(){
  return y_;
}

float NativeAtom::z(){
  return z_;
}

std::string NativeAtom::print() const{
  return "ATOM " + utility::to_string(type_) + " " + utility::to_string(x_) + " " + utility::to_string(y_) + " " +
      utility::to_string(z_) + "\n";
}

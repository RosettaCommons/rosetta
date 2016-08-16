// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /rosetta/rosetta_source/src/protocols/features/helixAssembly/HelicalFragment.cc
///
/// @brief

/// @author Tim jacobs

//Unit Headers
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Core
#include <core/conformation/Residue.hh>

//Utility
#include <utility/string_util.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>

//C++
#include <cmath>
#include <iostream>

namespace protocols {
namespace features {
namespace helixAssembly {

HelicalFragment::HelicalFragment():start_(0),end_(0){}

HelicalFragment::HelicalFragment(core::Size start, core::Size end):
	start_(start),
	end_(end)
{}

HelicalFragment::~HelicalFragment(){}

core::Size HelicalFragment::start() const { return start_; }
core::Size HelicalFragment::seq_start() const { return std::min(start_, end_); }

core::Size HelicalFragment::end() const { return end_; }
core::Size HelicalFragment::seq_end() const { return std::max(start_, end_); }
bool HelicalFragment::reversed() const { return start_ > end_; }

numeric::xyzVector<core::Real> HelicalFragment::com() const{ return com_; }
void HelicalFragment::com(numeric::xyzVector<core::Real> com)
{
	this->com_ = com;
}

numeric::xyzVector<core::Real> HelicalFragment::p0() const{ return p0_; }
void HelicalFragment::p0(numeric::xyzVector<core::Real> p0)
{
	this->p0_ = p0;
}

numeric::xyzVector<core::Real> HelicalFragment::p1() const{ return p1_; }
void HelicalFragment::p1(numeric::xyzVector<core::Real> p1)
{
	this->p1_ = p1;
}

numeric::xyzVector<core::Real> HelicalFragment::principal_component() const { return principal_component_; }
void HelicalFragment::principal_component(numeric::xyzVector<core::Real> principal_component)
{
	this->principal_component_ = principal_component;
}

core::Real HelicalFragment::sasa() const { return sasa_; }
void HelicalFragment::sasa(core::Real sasa)
{
	this->sasa_ = sasa;
}

core::Size HelicalFragment::size() const
{
	return (end_ > start_) ? end_-start_+1 : start_-end_+1;
}

std::ostream &
operator <<( std::ostream & os, HelicalFragment const & t )
{
	os << "(" << t.start() << ", " << t.end() << ")";
	return os;
}

} //namespace helixAssembly
} //namespace features
} //namespace protocols


// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionOptions.cc
/// @brief  Options for fiber diffraction data
/// @author Wojciech Potrzebowski and Ingemar Andre

// Unit Headers
#include <core/scoring/fiber_diffraction/FiberDiffractionOptionsCreator.hh>
#include <core/scoring/fiber_diffraction/FiberDiffractionOptions.hh>
#include <utility/tag/Tag.hh>

// Plaform Headers
#include <core/types.hh>

// C++ Headers
#include <string>

namespace core {
namespace scoring {
namespace fiber_diffraction {

using basic::resource_manager::ResourceOptionsOP;
using basic::resource_manager::ResourceOptions;
using std::string;
using utility::tag::TagCOP;

FiberDiffractionOptionsCreator::FiberDiffractionOptionsCreator() {}

FiberDiffractionOptionsCreator::~FiberDiffractionOptionsCreator() {}

ResourceOptionsOP
FiberDiffractionOptionsCreator::create_options() const {
	return ResourceOptionsOP( new FiberDiffractionOptions );
}

string
FiberDiffractionOptionsCreator::options_type() const {
	return "FiberDiffractionOptions";
}

//3A and 12A resolution cutoffs. 
//c repeat 3.0A is typical for inoviruses but has to be adjusted for each system
FiberDiffractionOptions::FiberDiffractionOptions() :
	ResourceOptions(),
	c(3.0),
	res_cutoff_high(0.333333),
	res_cutoff_low(0.0833333)	
{}


FiberDiffractionOptions::FiberDiffractionOptions(
	string const & name
) :
	ResourceOptions(name),
	c(3.0),
	res_cutoff_high(0.333333),
  res_cutoff_low(0.0833333)
{}

FiberDiffractionOptions::FiberDiffractionOptions(
  string const & name,
	Real c_,
	Real res_cutoff_high_,
	Real res_cutoff_low_
) :
  ResourceOptions(name),
	c(c_),
  res_cutoff_high(res_cutoff_high_),
  res_cutoff_low(res_cutoff_low_)
{}


FiberDiffractionOptions::FiberDiffractionOptions(
	FiberDiffractionOptions const & src
) :
	ResourceOptions(src),
	c(src.c),
  res_cutoff_high(src.res_cutoff_high),
  res_cutoff_low(src.res_cutoff_low)
{}

FiberDiffractionOptions::~FiberDiffractionOptions() {}

Real FiberDiffractionOptions::get_res_high() const {
	return  res_cutoff_high;
	}

void  FiberDiffractionOptions::set_res_high( Real res_cutoff_high_ ) {
	res_cutoff_high = res_cutoff_high_;
}

Real FiberDiffractionOptions::get_res_low() const {
  return  res_cutoff_low;
  }

void  FiberDiffractionOptions::set_res_low( Real res_cutoff_low_ ) {
  res_cutoff_low = res_cutoff_low_;
}

Real FiberDiffractionOptions::get_c_repeat() const {
  return c;
}

void  FiberDiffractionOptions::set_c_repeat( Real c_ ) {
  c = c_;
}

void
FiberDiffractionOptions::parse_my_tag(
  TagCOP tag
) {
	c = tag->getOption<Real>("c", 3.0);
	res_cutoff_high = tag->getOption<Real>("res_cutoff_high", 0.333333);
  res_cutoff_low = tag->getOption<Real>("res_cutoff_low", 0.0833333);

}


} // namespace
} // namespace
} // namespace

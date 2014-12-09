// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/SasaMethod.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/scoring/sasa/SasaMethod.hh>

namespace core {
namespace scoring {
namespace sasa {
	
SasaMethod::SasaMethod(Real probe_radius, SasaRadii radii_set):
	utility::pointer::ReferenceCount(),
	probe_radius_(probe_radius),
	radii_set_(radii_set),
	include_probe_radius_(true),
	use_big_polar_H_(false)
{
		
}

SasaMethod::~SasaMethod(){}

void
SasaMethod::set_include_probe_radius_in_calc(bool include_probe_radius) {
	include_probe_radius_ = include_probe_radius;
}

void
SasaMethod::set_probe_radius(Real probe_radius){
	probe_radius_ = probe_radius;
}

void
SasaMethod::set_radii_set(SasaRadii radii_set) {
	radii_set_ = radii_set;
}

void
SasaMethod::set_use_big_polar_hydrogen(bool big_polar_h){
	use_big_polar_H_ = big_polar_h;
}


}
} //scoring 
} //core

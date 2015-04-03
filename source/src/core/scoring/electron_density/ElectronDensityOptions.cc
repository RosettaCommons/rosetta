// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/electron_density/ElectronDensityOptions.cc
/// @brief  Options for constructing an electron density map
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <core/scoring/electron_density/ElectronDensityOptionsCreator.hh>
#include <core/scoring/electron_density/ElectronDensityOptions.hh>
#include <utility/tag/Tag.hh>

// Plaform Headers
#include <core/types.hh>

// C++ Headers
#include <string>

namespace core {
namespace scoring {
namespace electron_density {

using basic::resource_manager::ResourceOptionsOP;
using basic::resource_manager::ResourceOptions;
using std::string;
using utility::tag::TagCOP;

///// ElectronDensityOptionsCreator /////
ElectronDensityOptionsCreator::ElectronDensityOptionsCreator() {}

ElectronDensityOptionsCreator::~ElectronDensityOptionsCreator() {}

ResourceOptionsOP
ElectronDensityOptionsCreator::create_options() const {
	return ResourceOptionsOP( new ElectronDensityOptions );
}

string
ElectronDensityOptionsCreator::options_type() const {
	return "ElectronDensityOptions";
}

ElectronDensityOptions::ElectronDensityOptions() :
	ResourceOptions(),
	mapreso_(3.0),
	grid_spacing_(0.0) // if <= 0 then do not resize gride spacing
{}

ElectronDensityOptions::ElectronDensityOptions(
	string const & name
) :
	ResourceOptions(name),
	mapreso_(3.0),
	grid_spacing_(0.0) // if <= 0 then do not resize gride spacing
{}

ElectronDensityOptions::ElectronDensityOptions(
	string const & name,
	Real mapreso,
	Real grid_spacing
) :
	ResourceOptions(name),
	mapreso_(mapreso),
	grid_spacing_(grid_spacing) // if <= 0 then do not resize gride spacing
{}

ElectronDensityOptions::~ElectronDensityOptions() {}

ElectronDensityOptions::ElectronDensityOptions(
	ElectronDensityOptions const & src
) :
	ResourceOptions(src),
	mapreso_(src.mapreso_),
	grid_spacing_(src.grid_spacing_)
{}

Real
ElectronDensityOptions::get_mapreso(
) const {
	return mapreso_;
}

void
ElectronDensityOptions::set_mapreso(
	Real mapreso
) {
	mapreso_ = mapreso;
}

Real
ElectronDensityOptions::get_grid_spacing(
) const {
	return grid_spacing_;
}

void
ElectronDensityOptions::set_grid_spacing(
	Real grid_spacing
) {
	grid_spacing_ = grid_spacing;
}

void
ElectronDensityOptions::parse_my_tag(
	TagCOP tag
) {
	mapreso_ = tag->getOption<Real>("mapreso", 3.0);
	grid_spacing_ = tag->getOption<Real>("grid_spacing", 0.0);
}


} // namespace
} // namespace
} // namespace

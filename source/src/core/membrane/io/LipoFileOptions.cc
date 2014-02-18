// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileOptions.cc
///
/// @brief      Options Class - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipoFileOpitons_cc
#define INCLUDED_core_membrane_io_LipoFileOpitons_cc

// Unit Headers
#include <core/membrane/io/LipoFileOptions.hh>
#include <core/membrane/io/LipoFileOptionsCreator.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace membrane {
namespace io {

/// @brief Constructor
LipoFileOptions::LipoFileOptions() {}

/// @brief Destructor
LipoFileOptions::~LipoFileOptions() {}

/// @brief Read options from .xml file
void
LipoFileOptions::parse_my_tag(
	utility::tag::TagCOP
	)
{}

/// @brief Return options class type
std::string
LipoFileOptions::type() const
{
	return "LipoFileOptions";
}

/// @brief Creator class - return lipo file options class type
std::string
LipoFileOptionsCreator::options_type() const { return "LipoFileOptions"; }

/// @brief Creator class - return lipo file options class
basic::resource_manager::ResourceOptionsOP
LipoFileOptionsCreator::create_options() const { return new LipoFileOptions; }


} // io
} // membrane
} // core

#endif INCLUDED_core_membrane_io_LipoFileOpitons_cc

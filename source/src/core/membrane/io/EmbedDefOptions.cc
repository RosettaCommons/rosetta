// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefOptions.cc
///
/// @brief      Options for membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefOptions_cc
#define INCLUDED_core_membrane_io_EmbedDefOptions_cc

// Unit Headers
#include <core/membrane/io/EmbedDefOptions.hh>
#include <core/membrane/io/EmbedDefOptionsCreator.hh>

// Package Headers
#include <core/types.hh>

#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>
#include <string>

static basic::Tracer TR( "core.membrane.io.EmbedDefOptions" );

namespace core {
namespace membrane {
namespace io {

/// @brief Constructor
EmbedDefOptions::EmbedDefOptions() {}

/// @brief Destructor
EmbedDefOptions::~EmbedDefOptions() {}

/// @brief Parse .xml file for options
void
EmbedDefOptions::parse_my_tag(
	utility::tag::TagCOP
	)
{}

/// @brief Options class type embed def options
std::string
EmbedDefOptions::type() const
{
	return "EmbedDefOptions";
}

/// @brief Options type - embedding definition
std::string
EmbedDefOptionsCreator::options_type() const { return "EmbedDef"; }

/// @brief Return new options class embed def options
basic::resource_manager::ResourceOptionsOP
EmbedDefOptionsCreator::create_options() const { return new EmbedDefOptions; }

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefOptions_cc



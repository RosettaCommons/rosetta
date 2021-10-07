// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/bcl/util.hh
/// @brief Utilities for interacting with the BCL.
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Benjamin P. Brown (benjamin.p.brown17@gmail.com)

#ifndef INCLUDED_core_chemical_bcl_util_hh
#define INCLUDED_core_chemical_bcl_util_hh

#include <core/types.hh>

namespace core {
namespace chemical {
namespace bcl {

/// @brief Utility function to print BCL usage details and then utility_exit
void require_bcl();

/// @brief Initialize core components of the bcl apps/apps.cpp
/// @details e.g. GetExecutablePath, etc.
void initialize_bcl_main();

/// @brief Initialize the BCL random number generator.
/// @details Note that seed is an int to match the seed generated in core/init.cc
void initialize_bcl_random( int const seed );

/// @brief Initialize the BCL output levels with the Rosetta commandline settings
/// @details You can set the global BCL output by controlling the "BCL" tracer.
void initialize_bcl_tracers();

/// @brief Locate the BCL executable path within /main/source/externals/bcl
/// @details Required that we find the BCL submodule executable path so that
/// we can set paths to the rotamer library and model directories; this is based
/// obviously on the locate_database() function in init.cc with some extra logic
/// and minor additional differences
void locate_bcl();

} // namespace bcl
} // namespace chemical
} // namespace core

#endif // Include guard.

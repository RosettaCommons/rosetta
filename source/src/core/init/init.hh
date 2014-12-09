// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/init/init.hh
/// @brief  Core Init functions
/// @author Sergey Lyskov
///


#ifndef INCLUDED_core_init_init_hh
#define INCLUDED_core_init_init_hh

// STL
#include <string>

#include <utility/vector1.hh>
#include <numeric/random/random.hh>

namespace core {
namespace init {

/// @brief Initialize MPI (message passing interface) for parallel applications
void init_mpi(int argc, char * argv []);

/// @brief Initialize the option system, which manages command line options
void init_options(int argc, char * argv []);

/// @brief Tracers control output to std::cout and std::cerr
void init_tracers();

/// @brief Choose to output source version control information?
void init_source_revision();

/// @brief Setup basic search paths
void init_paths();

/// @brief Check for deprecated flags specified by the user and output error messages if necessary
void check_deprecated_flags();

/// @brief Describe the application execution command
void report_application_command(int argc, char * argv []);

/// @brief Initalize random number generators
void init_random_number_generators();

/// @brief Choose to randomly delay execution to desyncronize parallel execution
void random_delay();

/// @brief Locate rosetta_database
void locate_rosetta_database();

/// @brief Profiling measures execution performance
void init_profiling();

/// @brief Init basic core systems: options system, random system.
void init(int argc, char * argv []);

/// @brief wrapper for core system Init
void init( utility::vector1<std::string> const & args );

/// @brief Initialize random generator systems (and send debug io to tracer with seed/mode info).
void init_random_generators( int start_seed, std::string const & RGtype );

/// @brief Check for deprecated flags and utility exit if any deprecated flags are detected
void check_deprecated_flags();

} // namespace init
} // namespace core

#endif // INCLUDED_core_init_init_HH

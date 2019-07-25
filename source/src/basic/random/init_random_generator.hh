// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/random/init_random_generator.hh
/// @brief  Functions to initialize the random generator.
/// @author Originally by Sergey Lyskov.
/// @modified Moved to basic from core by Vikram K. Mulligan (vmulligan@flatironinstitute.org).  Added multithreading support
/// to ensure that all threads have unique random seeds.


#ifndef INCLUDED_basic_random_init_random_generator_hh
#define INCLUDED_basic_random_init_random_generator_hh

#include <string>
#include <basic/random/RandomGeneratorSettings.fwd.hh>
#include <platform/types.hh>

namespace basic {
namespace random {

/// @brief Initalize random number generators from the command line
/// @details If this is the multi-threaded build, threads other than thread zero must
/// provide their thread index to ensure that they get unique random seeds.
void init_random_number_generators(
#ifdef MULTI_THREADED
	platform::Size const curthread=0
#endif
);

/// @brief Figure out what seed to use based on a variety of settings stored in the
/// RandomGeneratorSettings object.
/// @note In MPI mode, this function uses MPI_Bcast calls, and that these calls
/// require that all MPI nodes reach this function, or deadlock will occur. This
/// deadlock can be avoided if the settings' mpi_bcast_seed_from_node0_ boolean is set
/// to false.
/// @details If this is the multi-threaded build, threads other than thread zero must
/// provide their thread index to ensure that they get unique random seeds.
int determine_random_number_seed(
	basic::random::RandomGeneratorSettings const & settings
#ifdef MULTI_THREADED
	,
	platform::Size const curthread
#endif
);

#ifdef MULTI_THREADED
/// @brief Given a seed that has already been adjusted for MPI (if this is the MPI build),
/// adjust for the current thread index.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void adjust_seed_for_thread( int & real_seed, platform::Size const curthread );
#endif

/// @brief Given a particular seed, initialize random generator systems (and send debug io to tracer with seed/mode info).
void init_random_generators( int start_seed, std::string const & RGtype );

} // namespace random
} // namespace basic

#endif // INCLUDED_basic_random_init_random_generator_hh

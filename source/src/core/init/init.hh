// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/init/init.hh
/// @brief  Core Init functions
/// @author Sergey Lyskov


#ifndef INCLUDED_core_init_init_hh
#define INCLUDED_core_init_init_hh

// STL
#include <string>

#include <utility/vector1.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
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

class RandomGeneratorSettings
{
public:
	RandomGeneratorSettings();
	void initialize_from_options( utility::options::OptionCollection const & options );
	static void list_options_read( utility::options::OptionKeyList & opt_keys );

	int seed() const;
	int seed_offset() const;
	bool const_seed() const;
	bool use_time_as_seed() const;
	std::string const & random_device_name() const;
	std::string const & rng_type() const;
	bool mpi_bcast_seed_from_node0() const;

	void seed( int setting );
	void seed_offset( int setting );
	void const_seed( bool setting );
	void use_time_as_seed( bool setting );
	void random_device_name( std::string const & setting );
	void rng_type( std::string const & setting );
	/// @brief It is important that calls to determine_random_seed
	/// with this value set to true (the default!) only occur when all
	/// nodes will reach this function at the same time with no intervening
	/// mpi calls. Beware if in a multithreaded context if a seed is set
	/// more than once (e.g. after each thread launches) as it can produce
	/// MPI-deadlock.
	void mpi_bcast_seed_from_node0( bool setting );

private:
	int seed_;
	int seed_offset_;
	bool const_seed_;
	bool use_time_as_seed_;
	std::string random_device_name_;
	std::string rng_type_;
	bool mpi_bcast_seed_from_node0_;
};

/// @brief Initalize random number generators from the command line
void init_random_number_generators();

/// @brief Figure out what seed to use based on a variety of settings stored in the
/// RandomGeneratorSettings object.
/// @note In MPI mode, this function uses MPI_Bcast calls, and that these calls
/// require that all MPI nodes reach this function, or deadlock will occur. This
/// deadlock can be avoided if the settings' mpi_bcast_seed_from_node0_ boolean is set
/// to false.
int determine_random_number_seed( RandomGeneratorSettings const & settings );

/// @brief Given a particular seed, initialize random generator systems (and send debug io to tracer with seed/mode info).
void init_random_generators( int start_seed, std::string const & RGtype );

/// @brief Choose to randomly delay execution to desyncronize parallel execution
void random_delay();

/// @brief Locate rosetta_database
void locate_rosetta_database();

/// @brief Profiling measures execution performance
void init_profiling();

/// @brief Set up system resources
void init_resources();

/// @brief Init basic core systems: options system, random system.
void init(int argc, char * argv []);

/// @brief wrapper for core system Init
void init( utility::vector1<std::string> const & args );

/// @brief Check for deprecated flags and utility exit if any deprecated flags are detected
void check_deprecated_flags();

} // namespace init
} // namespace core

#endif // INCLUDED_core_init_init_HH

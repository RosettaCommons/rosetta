// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/random/RandomGeneratorSettings.hh
/// @brief A class to store settings from the options system for the random generator.  Moved from core to basic.
/// @author Original author unknown.
/// @modified Moved from core to basic by Vikram K. Mulligan (vmulligan@flatironinstitute.org).

#ifndef INCLUDED_basic_random_RandomGeneratorSettings_hh
#define INCLUDED_basic_random_RandomGeneratorSettings_hh

#include <basic/random/RandomGeneratorSettings.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/options/OptionCollection.hh>

namespace basic {
namespace random {

/// @brief A class to store settings from the options system for the random generator.  Moved from core to basic.
/// @author Original author unknown.
/// @modified Moved from core to basic by Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class RandomGeneratorSettings : public utility::pointer::ReferenceCount
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

} //basic
} //random

#endif //INCLUDED_basic_random_RandomGeneratorSettings_hh

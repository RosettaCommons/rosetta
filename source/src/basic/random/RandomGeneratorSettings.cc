// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/random/RandomGeneratorSettings.cc
/// @brief A class to store settings from the options system for the random generator.  Moved from core to basic.
/// @author Original author unknown.
/// @modified Moved from core to basic by Vikram K. Mulligan (vmulligan@flatironinstitute.org).

// Project headers:
#include <basic/random/RandomGeneratorSettings.hh>

// Basic headers:
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <utility/options/keys/OptionKeyList.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "basic.random.RandomGeneratorSettings" );


namespace basic {
namespace random {

RandomGeneratorSettings::RandomGeneratorSettings() :
	seed_( 1111111 ),
	seed_offset_( 0 ),
	const_seed_( false ),
	use_time_as_seed_( false ),
	random_device_name_( "/dev/urandom" ),
	rng_type_( "mt19937" ),
	mpi_bcast_seed_from_node0_( true )
{}

void RandomGeneratorSettings::initialize_from_options( utility::options::OptionCollection const & options )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( options[ run::constant_seed ].active() )  const_seed_ = options[ run::constant_seed ]();
	if ( options[ run::jran ].active() )  seed_ = options[ run::jran ]();
	if ( options[ run::seed_offset ].active() )  seed_offset_ = options[ run::seed_offset ]();
	if ( options[ run::use_time_as_seed ].active() )  use_time_as_seed_ = options[ run::use_time_as_seed ]();

	random_device_name_ = options[ run::rng_seed_device ](); // typically /dev/urandom or /dev/random
	rng_type_ = options[ run::rng ];
}

void RandomGeneratorSettings::list_options_read( utility::options::OptionKeyList & opt_keys )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	opt_keys
		+ run::constant_seed
		+ run::jran
		+ run::seed_offset
		+ run::use_time_as_seed
		+ run::rng_seed_device
		+ run::rng;
}

int RandomGeneratorSettings::seed() const { return seed_; }
int RandomGeneratorSettings::seed_offset() const { return seed_offset_; }
bool RandomGeneratorSettings::const_seed() const { return const_seed_; }
bool RandomGeneratorSettings::use_time_as_seed() const { return use_time_as_seed_; }
std::string const & RandomGeneratorSettings::random_device_name() const { return random_device_name_; }
std::string const & RandomGeneratorSettings::rng_type() const { return rng_type_; }
bool RandomGeneratorSettings::mpi_bcast_seed_from_node0() const { return mpi_bcast_seed_from_node0_; }

void RandomGeneratorSettings::seed( int setting ) { seed_ = setting; }
void RandomGeneratorSettings::seed_offset( int setting ) { seed_offset_ = setting; }
void RandomGeneratorSettings::const_seed( bool setting ) { const_seed_ = setting; }
void RandomGeneratorSettings::use_time_as_seed( bool setting ) { use_time_as_seed_ = setting; }
void RandomGeneratorSettings::random_device_name( std::string const & setting ) { random_device_name_ = setting; }
void RandomGeneratorSettings::rng_type( std::string const & setting ) { rng_type_ = setting; }
void RandomGeneratorSettings::mpi_bcast_seed_from_node0( bool setting ) { mpi_bcast_seed_from_node0_ = setting; }

} //random
} //basic

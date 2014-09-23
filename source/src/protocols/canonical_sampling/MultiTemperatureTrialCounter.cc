// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Phil Bradley

// Unit Headers
#include <protocols/canonical_sampling/MultiTemperatureTrialCounter.hh>
#include <protocols/canonical_sampling/SimulatedTempering.hh>



// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh> //Pretty output.
#include <utility/io/ozstream.hh>

// Utility Headers
#include <basic/Tracer.hh>
static thread_local basic::Tracer tr( "protocols.canonical_sampling.MultiTemperatureTrialCounter" );

namespace protocols {
namespace canonical_sampling {

using namespace core;
using namespace ObjexxFCL;

MultiTemperatureTrialCounter::MultiTemperatureTrialCounter( TemperatureController const * temp_in ) :
	tempering_( temp_in )
{
	reset();
}

void
MultiTemperatureTrialCounter::set_temperature_observer( TemperatureController const * temp_in ) {
	tempering_ = temp_in;
	reset();
}


/////////////////////////////////////////////////////////////////////////////
void
MultiTemperatureTrialCounter::reset()
{
	runtime_assert( tempering_ != 0 );
	counters_.clear();
	counters_.resize( tempering_->n_temp_levels() );
}

void MultiTemperatureTrialCounter::show() const {
	show( tr.Info );
}
/////////////////////////////////////////////////////////////////////////////
void
MultiTemperatureTrialCounter::show( std::ostream& os ) const {
	assert( tempering_ );
	_write_to_stream( os, "" );
}

void
MultiTemperatureTrialCounter::write_to_file( std::string const& file, std::string const& tag ) const {
	utility::io::ozstream out( file, std::ios::app );
	_write_to_stream( out, tag );
}

void
MultiTemperatureTrialCounter::_write_to_stream( std::ostream& os, std::string const& tag ) const {
	assert( tempering_ );
	os << "Stats for job: " << tag << std::endl;
	for ( Size i=1; i<=counters_.size(); ++i ) {
		std::ostringstream line_header;
		line_header << "level " << format::I( 4, i) << " temperature " << format::F( 4, 2, tempering_->temperature( i ) );
		counters_[ i ].show( os, line_header.str(), true );
	}
}

void MultiTemperatureTrialCounter::count_trial( std::string const& tag ) {
	assert( tempering_ );
	assert( tempering_->n_temp_levels() == counters_.size() );
	counters_[ tempering_->temperature_level() ].count_trial( tag );
}

void MultiTemperatureTrialCounter::count_accepted( std::string const& tag ) {
	assert( tempering_ );
	counters_[ tempering_->temperature_level() ].count_accepted( tag );
}

void MultiTemperatureTrialCounter::count_energy_drop( std::string const& tag, Real delta ) {
	assert( tempering_ );
	counters_[ tempering_->temperature_level() ].count_energy_drop( tag, delta );
}

protocols::moves::TrialCounter const&
MultiTemperatureTrialCounter::operator[]( core::Size level ) const {
	return counters_[ level ];
}

protocols::moves::TrialCounter&
MultiTemperatureTrialCounter::operator[]( core::Size level ) {
	return counters_[ level ];
}

} // namespace canonical_sampling
} // namespace core

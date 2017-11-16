// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Phil Bradley

// Unit Headers
#include <protocols/moves/TrialCounter.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh> //Pretty output.

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>

static basic::Tracer tr( "protocols.moves.TrialCounter" );

namespace protocols {
namespace moves {

using namespace core;
using namespace ObjexxFCL::format;

/////////////////////////////////////////////////////////////////////////////
void
TrialCounter::reset()
{
	trial_counter_.clear();
	accept_counter_.clear();
	energy_drop_counter_.clear();
}

/// @detail return number of trials since last reset
core::Size
TrialCounter::total_trials() const {
	Size ntrials( 0 );
	for ( auto const & it : trial_counter_ ) {
		ntrials += it.second;
	}
	return ntrials;
}

void TrialCounter::show() const {
	show( tr.Info );
}
/////////////////////////////////////////////////////////////////////////////
void
TrialCounter::show( std::ostream& os, std::string line_header, bool endline ) const {
	if ( line_header.size() ) {
		line_header=line_header+" ";
	}
	for ( auto const & it : trial_counter_ ) {
		std::string const & tag( it.first );
		int const ntrials( it.second );
		if ( accept_counter_.count( tag ) ) {
			int const accepts( accept_counter_.find( tag )->second );
			os << line_header << A( 16, tag ) <<
				" trials= " << I( 6, ntrials ) << "; " <<
				" accepts= " << F( 6, 4, core::Real( accepts )/ntrials ) << "; ";
			auto edc_it = energy_drop_counter_.find( tag );
			if ( edc_it != energy_drop_counter_.end() ) {
				core::Real const energy_drop( edc_it->second );
				os << " energy_drop/trial= " << F( 9, 5, core::Real( energy_drop ) / ntrials );
			}
			if ( endline ) os << std::endl;
			else os << " ";
		} else { //no accepts
			os << line_header << A( 16, tag ) << " trials= " << I( 6, ntrials ) << " NO ACCEPTS.";
			if ( endline ) os << std::endl;
			else os << " ";
		} // else
	} // for
}

utility::vector1< std::string > const TrialCounter::tags() const {
	utility::vector1< std::string > tags;
	std::map< std::string, int >::const_iterator it;
	for ( it = trial_counter_.begin(); it != trial_counter_.end(); ++it ) {
		tags.push_back(it->first);
	}
	return tags;
}

core::Size TrialCounter::trial( std::string const& tag ) const {
	std::map< std::string, int >::const_iterator it;
	it = trial_counter_.find( tag );
	return (it != trial_counter_.end()) ? it->second : 0;
}

core::Size TrialCounter::accepted( std::string const& tag ) const {
	std::map< std::string, int >::const_iterator it;
	it = accept_counter_.find( tag );
	return (it != accept_counter_.end()) ? it->second : 0;
}

core::Real TrialCounter::energy_drop( std::string const& tag ) const {
	std::map< std::string, core::Real >::const_iterator it;
	it = energy_drop_counter_.find( tag );
	return (it != energy_drop_counter_.end()) ? it->second : 0;
}

void TrialCounter::count_trial( std::string const& tag ) {
	++trial_counter_[ tag ];
}

void TrialCounter::count_accepted( std::string const& tag ) {
	++accept_counter_[ tag ];
}

void TrialCounter::count_energy_drop( std::string const& tag, Real delta ) {
	energy_drop_counter_[ tag ] += delta;
}


} // namespace moves
} // namespace core

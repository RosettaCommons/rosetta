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
#include <protocols/moves/TrialCounter.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh> //Pretty output.

// Utility Headers
#include <basic/Tracer.hh>
static basic::Tracer tr("protocols.moves.TrialCounter");

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

///@detail return number of trials since last reset
core::Size
TrialCounter::total_trials() const {
	Size ntrials( 0 );
	for ( std::map< std::string, int >::const_iterator
					it=trial_counter_.begin(); it != trial_counter_.end(); ++it ) {
		ntrials += it->second;
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
	for ( std::map< std::string, int >::const_iterator
					it=trial_counter_.begin(); it != trial_counter_.end(); ++it ) {
		std::string const & tag( it->first );
		int const ntrials( it->second );
		if ( accept_counter_.count( tag )) {
			int const accepts( accept_counter_.find( tag )->second );
			core::Real const energy_drop( energy_drop_counter_.find( tag )->second );
			os << line_header << A( 16, tag ) <<
				" trials= " << I( 6, ntrials ) << "; " <<
				" accepts= " << F( 6, 4, core::Real( accepts )/ntrials ) << "; " <<
				" energy_drop/trial= " << F( 9, 5, core::Real( energy_drop ) / ntrials );
			if ( endline ) os << std::endl;
		} else { //no accepts
			os << A( 16, tag ) << " trials= " << I( 6, ntrials ) <<	" NO ACCEPTS.";
			if ( endline ) os << std::endl;
		} // else
	} // for
}

core::Size TrialCounter::trial( std::string const& tag ) {
	return trial_counter_[ tag ];
}

core::Size TrialCounter::accepted( std::string const& tag ) {
	return accept_counter_[ tag ];
}

core::Real TrialCounter::energy_drop( std::string const& tag ) {
	return energy_drop_counter_[ tag ];
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

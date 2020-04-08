// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/util.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/moves/util.hh>

// Utility headers
#include <utility/exit.hh>

// Project headers
#include <protocols/filters/Filter.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {

protocols::moves::MoverOP find_mover_or_die(
	const std::string& mover_name,
	const utility::tag::TagCOP,
	const basic::datacache::DataMap & data
) {
	if ( ! data.has( "movers", mover_name ) ) {
		utility_exit_with_message(mover_name + " was not found in the list of availible Movers");
	}
	return data.get_ptr< protocols::moves::Mover >( "movers", mover_name );
}

protocols::filters::FilterOP find_filter_or_die(
	const std::string& filter_name,
	const utility::tag::TagCOP,
	const basic::datacache::DataMap & data
) {
	if ( ! data.has( "filters", filter_name ) ) {
		utility_exit_with_message(filter_name + " was not found in the list of availible Filters");
	}
	return data.get_ptr< protocols::filters::Filter >( "filters", filter_name );
}


}  // namespace moves
}  // namespace protocols

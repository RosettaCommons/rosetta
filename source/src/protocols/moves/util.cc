// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {

protocols::moves::MoverOP find_mover_or_die(const std::string& mover_name,
	const utility::tag::TagCOP,
	const protocols::moves::Movers_map& movers) {
	protocols::moves::Movers_map::const_iterator i = movers.find(mover_name);
	if ( i == movers.end() ) {
		utility_exit_with_message(mover_name + " not found in Movers_map");
	}

	return i->second;
}

protocols::filters::FilterOP find_filter_or_die(const std::string& filter_name,
	const utility::tag::TagCOP,
	const protocols::filters::Filters_map& filters) {
	protocols::filters::Filters_map::const_iterator i = filters.find(filter_name);
	if ( i == filters.end() ) {
		utility_exit_with_message(filter_name + " not found in Filters_map");
	}

	return i->second;
}

}  // namespace moves
}  // namespace protocols

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/util.hh
/// @brief Util functions for filters.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_filters_util_hh
#define INCLUDED_protocols_filters_util_hh

#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/BasicFilters.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <vector>
#include <map>
#include <string>

namespace protocols {
namespace filters {

std::vector< std::pair< FilterOP, boolean_operations > >
create_compound_statement_from_logic(std::string const & name, basic::datacache::DataMap const & data);

///@brief Parses a name and splits it into a basic CompoundFilter OR returns the found filter.
FilterOP
parse_filter_with_logic(std::string const & name, basic::datacache::DataMap const & data);

static const std::map< std::string, filters::boolean_operations >
boolean_types(
	{
	{ "AND", AND  },
	{ "OR", OR },
	{ "XOR", XOR },
	{ "NOR", NOR },
	{ "NAND", NAND },
	{ "ORNOT", ORNOT },
	{ "ANDNOT", ANDNOT },
	{ "NOT", NOT }
	}
);

} //filters
} //protocols


#endif //protocols/filters_util_hh


// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/util.cc
/// @brief Util functions for filters.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/filters/util.hh>
#include <protocols/filters/BasicFilters.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

static basic::Tracer TR( "protocols.filters.util" );


namespace protocols {
namespace filters {

std::vector< std::pair< FilterOP, boolean_operations > >
create_compound_statement_from_logic(std::string const & name, basic::datacache::DataMap const & data){
	utility::vector1< std::string > nameSP = utility::string_split(name, ' ');
	std::vector< std::pair< FilterOP, boolean_operations > > statement;
	//Parse the logic here.
	//Can start with filter and filter; which we will need to include as and.
	utility::vector1< std::string > full_vec;
	if ( ! boolean_types.count(utility::uppercased(nameSP[1])) ) {
		if ( ! data.has( "filters", nameSP[1] ) ) {
			utility_exit_with_message("filter name or boolean operation not recognized: "+nameSP[1]);
		} else {
			full_vec.push_back("AND");
		}

	}
	for ( auto s : nameSP ) {
		full_vec.push_back(s);
	}
	if ( full_vec.size() % 2 != 0 ) {
		TR << "Passed Logic String: " << name << std::endl;
		utility_exit_with_message("Input logic for filter parsing is missing a field! ");
	}

	for ( core::Size i = 1; i <= core::Size(full_vec.size()/2); ++i ) {
		core::Size bool_position;
		core::Size filter_position;
		if ( i == 1 ) {
			bool_position = 1;
			filter_position = 2;
		} else {
			bool_position = (i * 2) - 1;
			filter_position = i * 2;
		}
		TR << "Parsing: " << full_vec[bool_position] <<" " << full_vec[filter_position] << std::endl;

		std::string boolean_val = utility::uppercased(full_vec[ bool_position]);
		if ( ! boolean_types.count(boolean_val) ) {
			utility_exit_with_message("Boolean value not understood! "+boolean_val);
		}
		boolean_operations bool_op = boolean_types.at(boolean_val);
		FilterOP filter = data.get_ptr< protocols::filters::Filter >( "filters", full_vec[filter_position] );
		if ( filter == nullptr ) {
			utility_exit_with_message("Filter not found in datamap! "+full_vec[filter_position]);
		}

		std::pair< FilterOP, boolean_operations > filt_pair = std::make_pair( filter, bool_op);
		statement.push_back( filt_pair);
	}
	return statement;
}

FilterOP
parse_filter_with_logic(std::string const & name, basic::datacache::DataMap const & data){
	utility::vector1< std::string > nameSP = utility::string_split(name, ' ');
	if ( nameSP.size() == 0 ) {
		utility_exit_with_message("No filter given in name!");
	}

	if ( nameSP.size() == 1 ) {
		std::string filter_name = nameSP[1];
		if ( ! data.has( "filters", filter_name ) ) {
			utility_exit_with_message(filter_name + " was not found in the list of availible Filters");
		}
		return data.get_ptr< protocols::filters::Filter >( "filters", filter_name );
	} else {
		std::vector< std::pair< FilterOP, boolean_operations > > statement = create_compound_statement_from_logic(name, data);
		CompoundFilterOP cp_filter = utility::pointer::make_shared< CompoundFilter >( statement);
		return cp_filter;
	}
	return nullptr;
}




} //filters
} //protocols



// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DatabaseEntryWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

//Unit
#include <protocols/wum/DatabaseEntryWorkUnit.hh>

//Basic
#include <basic/Tracer.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

//C++
#include <string>
#include <map>

static basic::Tracer TR( "DatabaseEntryWorkUnit" );

namespace protocols {
namespace wum {

using namespace std;

DatabaseEntryWorkUnit::DatabaseEntryWorkUnit(utility::sql_database::sessionOP db_session):
	db_session_(std::move(db_session))
{}

DatabaseEntryWorkUnit::DatabaseEntryWorkUnit( std::map<std::string,std::string> const & row_map ):
	WorkUnitBase(),
	row_map_(row_map)
{}

void
DatabaseEntryWorkUnit::serialize(){

	TR << "Serializing db entry data" << endl;

	//serialize the row map using commas to separate column name and value and a forward slash to separate columns
	string data("");
	for ( map<string,string>::const_iterator it = row_map_.begin();
			it != row_map_.end(); ++it ) {

		data += it->first + "," + it->second + "/";
	}

	//Add the results_query_string_ data to the end of the serial data separate by a pipe
	data+="|" + result_query_string_;
	serial_data() = data;
}

void
DatabaseEntryWorkUnit::deserialize(){
	const std::string data(serial_data());

	TR << "De-serializing db entry data" << endl;

	//split between the serialized row map data and the query string
	utility::vector1< std::string > tokens = utility::string_split(data, '|');
	if ( tokens.size() != 2 ) {
		utility_exit_with_message("Error: DatabaseEntryWorkUnit failed to deserialize");
	}
	result_query_string_=tokens[2];

	//deserialize the map
	utility::vector1< std::string > entries = utility::string_split(tokens[1], '/');

	TR << "Total columns: " << entries.size() << endl;

	for ( Size i=0; i<entries.size(); ++i ) {
		if ( !entries[i].empty() ) {
			utility::vector1< std::string > key_values = utility::string_split(entries[i], ',');
			if ( key_values.size() != 2 ) {
				utility_exit_with_message("Error: DatabaseEntryWorkUnit failed to deserialize the results map");
			}
			row_map_[key_values[0]] = key_values[1];
		}
	}
}

}//namespace wum
}//namespace protocols

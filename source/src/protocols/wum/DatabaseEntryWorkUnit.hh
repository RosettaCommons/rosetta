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
/// @brief A work unit base class that can serialize and deserialize a map representing a database row (keys are column names, values are column values)
///        and a string representing a query to be executed by the master node upon completion of the workunit. This work unit should be treated as a
///        pure virtual since no run() function is implemented.

/// @author Tim Jacobs

#ifndef INCLUDED_protocols_wum_DatabaseEntryWorkUnit_hh
#define INCLUDED_protocols_wum_DatabaseEntryWorkUnit_hh

//Unit
#include <protocols/wum/DatabaseEntryWorkUnit.fwd.hh>
#include <protocols/wum/WorkUnitBase.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//C++
#include <string>
#include <map>

namespace protocols {
namespace wum {

class DatabaseEntryWorkUnit : public protocols::wum::WorkUnitBase {
public:

	DatabaseEntryWorkUnit(utility::sql_database::sessionOP db_session);

	DatabaseEntryWorkUnit( std::map<std::string,std::string> row_map );

	~DatabaseEntryWorkUnit() override= default;

	protocols::wum::WorkUnitBaseOP clone() const override {
		return protocols::wum::WorkUnitBaseOP( new DatabaseEntryWorkUnit( *this ) );
	}

	/// @brief Accessor for database query string
	std::string result_query_string(){return result_query_string_;}

protected:
	//    void set_defaults();

	/// @brief Serialize the row_map_
	void serialize() override;

	/// @brief Deserialize the row_map_
	void deserialize() override;

protected:

	/// @brief The database connection
	utility::sql_database::sessionOP db_session_;

	/// @brief map that represents a database row - keys are columns, values are values
	std::map<std::string,std::string> row_map_;

	/// @brief A string that stores the database query you want to run when finished with the work unit
	std::string result_query_string_;

};

}//namespace wum
}//namespace protocols

#endif

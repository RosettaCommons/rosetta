// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/message_listening/DatabaseSchemaGeneratorListener.cc
///
/// @brief manage the one-time initializiation of tables in a database
/// @author Matthew O'Meara (mattjomeara@gmail.com)



#include <basic/message_listening/MessageListener.fwd.hh>
#include <basic/message_listening/DatabaseSchemaGeneratorListener.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/pointer/ReferenceCount.hh>


#include <string>

namespace basic {
namespace message_listening {

using std::string;
using utility::sql_database::sessionOP;

std::string const TABLE_EXISTS("TABLE_EXISTS");
std::string const TABLE_DOES_NOT_EXIST("TABLE_DOES_NOT_EXIST");

void
DatabaseSchemaGeneratorListener::receive(
	string const & table_name
) {
	table_names_.insert(table_name);
}


///@details
/// Answer the question, "does the table with the given name exist?"
/// Input: table_name
/// Output: There is a table with that name in the database
bool
DatabaseSchemaGeneratorListener::request(
	std::string const & table_name,
	std::string & return_data
) {

	if(table_names_.find(table_name) != table_names_.end()){
		return_data = TABLE_EXISTS;
		return false;
	}

	return_data = TABLE_DOES_NOT_EXIST;
	return true;
}




} //namespace
} //namespace

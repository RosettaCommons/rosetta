// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 utility/sql_database/types.cc
/// @author Matthew O'Meara

#include <utility/sql_database/types.hh>
#include <utility/exit.hh>
#include <platform/types.hh>
#include <string>
#include <sstream>

namespace utility {
namespace sql_database {

using std::string;
using std::stringstream;

DatabaseMode::e
database_mode_from_name(
	std::string database_mode
) {
	if(!database_mode.compare("sqlite3")){
		return DatabaseMode::sqlite3;
	} else if(!database_mode.compare("mysql")){
		return DatabaseMode::mysql;
	} else if(!database_mode.compare("postgres")){
		return DatabaseMode::postgres;
	} else {
		utility_exit_with_message(
			"Unrecognized database mode: '" + database_mode + "'");
	}
	return DatabaseMode::sqlite3; // make compiler happy
}

TransactionMode::e
transaction_mode_from_name(
	std::string transaction_mode
){
	if(!transaction_mode.compare("none")){
		return TransactionMode::none;
	} else if(!transaction_mode.compare("standard")){
		return TransactionMode::standard;
	} else if(!transaction_mode.compare("chunk")){
		return TransactionMode::chunk;
	} else {
		utility_exit_with_message(
			"Unrecognized transaction mode: '" + transaction_mode + "'");
	}
	return TransactionMode::standard; // make compiler happy
}

std::string
name_from_transaction_mode(
	TransactionMode::e transaction_mode
){
	switch(transaction_mode){
	case TransactionMode::none:
		return "none";
	case TransactionMode::standard:
		return "standard";
	case TransactionMode::chunk:
		return "chunk";
	default:
		stringstream err_msg;
		err_msg
			<< "Unrecognized transaction mode: '"
			<< static_cast<platform::Size>(transaction_mode) << "'";
		utility_exit_with_message(err_msg.str());
#ifdef WIN32
		return "Error";
#endif
	}
}
	
std::string
name_from_database_mode(
	DatabaseMode::e database_mode
){
	switch(database_mode){
	case DatabaseMode::sqlite3:
		return "sqlite3";
	case DatabaseMode::mysql:
		return "mysql";
	case DatabaseMode::postgres:
		return "postgres";
	default:
		stringstream err_msg;
		err_msg
			<< "Unrecognized databse mode: '"
			<< static_cast<platform::Size>(database_mode) << "'"; 
		utility_exit_with_message(err_msg.str());
#ifdef WIN32
		return "Error";
#endif
	}
}

}
}

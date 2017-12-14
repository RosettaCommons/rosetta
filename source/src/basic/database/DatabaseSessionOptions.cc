// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/database/DatabaseSessionOptions.cc
/// @brief  load the database session
/// @author Tim Jacobs

//unit headers
#include <basic/database/DatabaseSessionOptions.hh>
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/database/sql_utils.hh>

//utility headers
#include <utility/tag/Tag.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//C++ headers
#include <istream>

namespace basic {
namespace database {

DatabaseSessionOptions::DatabaseSessionOptions() = default;

DatabaseSessionOptions::DatabaseSessionOptions(
	std::string const & name) :
	ResourceOptions(name)
{}

DatabaseSessionOptions::~DatabaseSessionOptions() = default;

void
DatabaseSessionOptions::parse_my_tag(
	utility::tag::TagCOP tag
){
	database_session_ = parse_database_connection(tag);
}

utility::sql_database::sessionOP
DatabaseSessionOptions::database_session() const { return database_session_; }

std::string
DatabaseSessionOptions::type() const { return "DatabaseSessionOptions"; }

} // database
} // basic

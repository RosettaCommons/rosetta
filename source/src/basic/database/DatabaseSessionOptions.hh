// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/database/DatabaseSessionOptions.hh
/// @brief  load the database session
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_DatabaseSessionOptions_hh
#define INCLUDED_basic_database_DatabaseSessionOptions_hh

//unit headers
#include <basic/resource_manager/ResourceOptions.hh>

//utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

//C++ headers
#include <istream>

namespace basic {
namespace database {

class DatabaseSessionOptions : public basic::resource_manager::ResourceOptions {
public:

	DatabaseSessionOptions();

	DatabaseSessionOptions(std::string const & name);

	virtual ~DatabaseSessionOptions();

	virtual
	void
	parse_my_tag(utility::tag::TagCOP);

	utility::sql_database::sessionOP
	database_session() const;

	std::string
	type() const;

private:

	utility::sql_database::sessionOP database_session_;
};

} // database
} // basic

#endif // include guard

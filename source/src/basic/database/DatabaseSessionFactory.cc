// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  load the database session
/// @file   basic/database/DatabaseSessionFactory.cc
/// @author Tim Jacobs

// Unit headers
#include <basic/database/DatabaseSessionFactory.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/tracer.hh>


namespace basic {
namespace database {

static thread_local basic::Tracer TR( "basic.database.databaseSessionFactory" );

utility::sql_database::sessionOP
DatabaseSessionFactory::create_database_session()
{
	return new utility::sql_database::session();
}

} // database
} // basic

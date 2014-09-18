// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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

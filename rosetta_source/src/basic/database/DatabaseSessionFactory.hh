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
/// @file   basic/database/DatabaseSessionFactor.hh
/// @author Tim Jacobs


#ifndef INCLUDED_basic_database_DatabaseSessionFactory_hh
#define INCLUDED_basic_database_DatabaseSessionFactory_hh

// Unit headers
#include <utility/sql_database/DatabaseSessionManager.hh>

namespace basic {
namespace database {

class DatabaseSessionFactory
{
	public:

	utility::sql_database::sessionOP
	create_database_session();
};

} // database
} // basic
#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/sql_database/DatabaseConnectionManager.fwd.hh
/// @brief  utility::sql_database::DatabaseConnectionManager forward declarations
/// @author Matthew O'Meara (mattjomear@gmail.com)
/// @author Sam Deluca
/// @author Chris Miles

#ifndef INCLUDED_utility_sql_database_DatabaseSessionManager_fwd_hh
#define INCLUDED_utility_sql_database_DatabaseSessionManager_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace sql_database {


class session;
typedef pointer::shared_ptr< session > sessionOP;
typedef pointer::shared_ptr< session const > sessionCOP;

class transaction;
typedef pointer::shared_ptr< transaction > transactionOP;
typedef pointer::shared_ptr< transaction const > transactionCOP;

class DatabaseConnectionManager;

}
}


#endif

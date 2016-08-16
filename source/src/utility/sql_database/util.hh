// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/sql_database/util.hh
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_utility_sql_database_util_hh
#define INCLUDED_utility_sql_database_util_hh

#ifdef USEMYSQL // for variables that are only used with mysql
	#define MYSQL_ONLY(x) x
#else
#define MYSQL_ONLY(x)
#endif


#ifdef USEPOSTGRES // for variables that are only used with postgres
	#define POSTGRES_ONLY(x) x
#else
#define POSTGRES_ONLY(x)
#endif

#endif // INCLUDED_utility_sql_database_util_hh


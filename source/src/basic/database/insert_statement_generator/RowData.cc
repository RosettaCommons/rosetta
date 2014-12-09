// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/insert_statement_generator/RowData.cc
///
/// @brief Column data object for the insert statement generator
/// @author Sam DeLuca


#include <basic/database/insert_statement_generator/RowData.hh>
#include <utility/string_util.hh>
#include <boost/uuid/uuid_io.hpp>
#include <utility/sql_database/DatabaseSessionManager.hh>

namespace basic {
namespace database {
namespace insert_statement_generator {

RowDataBase::RowDataBase(std::string const & column_name) : column_name_(column_name)
{

}

RowDataBase::~RowDataBase() {}

std::string RowDataBase::get_column_name() const
{
	return column_name_;
}



}
}
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Constraints.cc
///
/// @brief Constraints class for the schema generator framework
/// @author Tim Jacobs

#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/Column.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

#include <platform/types.hh>

#include <string>
#include <sstream>
#include <utility/exit.hh>

namespace basic {
namespace database {
namespace schema_generator {

/// @details Auto-generated virtual destructor
Constraint::~Constraint() {}

using std::string;
using std::stringstream;
using platform::Size;
using platform::Real;


Constraint::Constraint(
	Column column
) :
	columns_()
{
	columns_.push_back(column);
}

Constraint::Constraint(
	Columns columns
) :
	columns_(columns)
{}


UniqueConstraint::UniqueConstraint(
	Column column
) :
	Constraint(column)
{}

UniqueConstraint::UniqueConstraint(
	Columns columns
) :
	Constraint(columns)
{}

std::string
UniqueConstraint::print(
	utility::sql_database::sessionOP /*database_session*/
) const {
	stringstream constraint;

	constraint << "UNIQUE (";

	for ( Size i=1; i<=columns_.size(); ++i ) {
		if ( i!=1 ) {
			constraint << ", ";
		}
		constraint << columns_[i].name();
	}
	constraint << ")";
	return constraint.str();
}


GreaterThanConstraint::GreaterThanConstraint(
	Column column,
	platform::Real value
) :
	Constraint(column),
	value_(value)
{}

std::string
GreaterThanConstraint::print(
	utility::sql_database::sessionOP /*database_session*/
) const {
	stringstream constraint;

	assert(columns_.size() == 1);

	constraint << "CONSTRAINT " << columns_[1].name() << "_greater_than CHECK ";
	constraint << "(" << columns_[1].name() << " >= " << value_ << ")";
	return constraint.str();
}


} // schema_generator
} // namespace database
} // namespace utility


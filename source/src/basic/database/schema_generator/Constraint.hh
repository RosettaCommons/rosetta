// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/Constraints.hh
///
/// @brief Constraint class for the schema generator framework
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_schema_generator_Constraint_HH
#define INCLUDED_basic_database_schema_generator_Constraint_HH

#include <utility/pointer/ReferenceCount.hh>
#include <basic/database/schema_generator/Constraint.fwd.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

#include <basic/database/schema_generator/Column.hh>

#include <string>
#include <platform/types.hh>

namespace basic{
namespace database{
namespace schema_generator{

class Constraint : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Constraint();

	Constraint(Column column);

	Constraint(Columns columns);

	virtual std::string
	print(
		utility::sql_database::sessionOP
	) const=0;

protected:

	Columns columns_;
};

//class ComparisonConstraint : public Constraint
//{
//public:
//
//	enum ComparisonOperator{
//		GREATER_THAN,
//		LESS_THAN,
//		EQUAL_TO,
//		NOT_EQUAL_TO
//	};
//
//	ComparisonConstraint(Column column1, int comparison_operator, Column column2);
//
//	ComparisonConstraint(Column column, int comparator, std::string value);
//
//	ComparisonConstraint(Column, Column);
//
//	std::string print();
//
//private:
//	std::string string_value;
//
//};

class UniqueConstraint : public Constraint {
public:

	UniqueConstraint(Column column);
	UniqueConstraint(Columns columns);
	
	std::string
	print(
		utility::sql_database::sessionOP
	) const;
};

class GreaterThanConstraint : public Constraint {
public:

	GreaterThanConstraint(
		Column column,
		platform::Real value);

	std::string
	print(
		utility::sql_database::sessionOP
	) const;
	
private:
	platform::Real value_;
};

} // schema_generator
} // namespace database
} // namespace utility
#endif

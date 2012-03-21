// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Constraints.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_utility_sql_database_Constraint_HH
#define INCLUDED_utility_sql_database_Constraint_HH

#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/Constraint.fwd.hh>
#include <utility/vector1.hh>

#include <utility/sql_database/Column.hh>

#include <string>

namespace utility{
namespace sql_database{        

class Constraint : public utility::pointer::ReferenceCount  
{
public:
    
    Constraint(Column column);
    
    Constraint(utility::vector1<Column> columns);
    
    virtual std::string print() = 0;
    
protected:
    
    utility::vector1<Column> columns_;
};

//class ComparisonConstraint : public Constraint    
//{
//public:
//    
//    enum ComparisonOperator{
//        GREATER_THAN,
//        LESS_THAN,
//        EQUAL_TO,
//        NOT_EQUAL_TO
//    };
//    
//    ComparisonConstraint(Column column1, int comparison_operator, Column column2);
//    
//    ComparisonConstraint(Column column, int comparator, std::string value);
//    
//    ComparisonConstraint(Column, Column);
//    
//    std::string print();
//    
//private:
//    std::string string_value;
//    
//};

class UniqueConstraint : public Constraint {
public:
          
    UniqueConstraint(Column column);    
    UniqueConstraint(utility::vector1<Column> columns);    
    virtual std::string print();
    
};

    
} // namespace sql_database
} // namespace utility
#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Constraints.cc
///
/// @brief
/// @author Tim Jacobs

#include <utility/sql_database/Constraint.hh>
#include <utility/sql_database/Column.hh>

#include <string>
#include <utility/exit.hh>

namespace utility{
namespace sql_database{            
    
Constraint::Constraint(Column column){
    this->columns_.push_back(column);
}

Constraint::Constraint(utility::vector1<Column> columns):
columns_(columns)
{}
        
//    std::string ComparisonConstraint::print(){
//        std::string constraint_string;
//        
//        
//    }

UniqueConstraint::UniqueConstraint(Column column):
Constraint(column)
{}

UniqueConstraint::UniqueConstraint(utility::vector1<Column> columns):
Constraint(columns)
{}
        
std::string UniqueConstraint::print(){
    std::string constraint_string = "UNIQUE (";
    
    for(int i=1; i<=columns_.size(); ++i){
        if(i!=1){
            constraint_string += ",";
        }
        constraint_string += columns_[i].name();
    }
    constraint_string += ")";
    return constraint_string;
}
    
} // namespace sql_database
} // namespace utility
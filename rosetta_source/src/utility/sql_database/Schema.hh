// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Schema.hh
///
/// @brief
/// @author tim



#ifndef INCLUDED_utility_sql_database_Schema_HH
#define INCLUDED_utility_sql_database_Schema_HH

#include <utility/sql_database/PrimaryKey.hh>
#include <utility/sql_database/ForeignKey.hh>
#include <utility/sql_database/Column.hh>
#include <utility/sql_database/Constraint.hh>

#include <boost/unordered_set.hpp>

#include <string>
#include <set>

namespace utility{
namespace sql_database{

class Schema
{

public:
    
    Schema(std::string table_name);
    
    Schema(std::string table_name, PrimaryKey primary_key);
       
    void init();
    
    void add_foreign_key(ForeignKey key);
    
    void add_column(Column column);
    
    void add_constraint(ConstraintOP constraint);
    
    std::string print();
    
private:
    
    std::string database_mode_;
    
    std::string table_name_;
    PrimaryKey primary_key_;
    utility::vector1<Column> columns_;
    utility::vector1<ForeignKey> foreign_keys_;
    utility::vector1<ConstraintOP> constraints_;
};

} // namespace sql_database
} // namespace utility
    
#endif

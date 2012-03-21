// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PrimaryKey.cc
///
/// @brief
/// @author tim

#include <utility/sql_database/PrimaryKey.hh>
//#include <utility/sql_database/Column.hh>

//C++ Headers
#include <string>
#include <set>

namespace utility{
namespace sql_database{

PrimaryKey::PrimaryKey(){}
    
PrimaryKey::PrimaryKey(Column column){
    columns_.push_back(column);
}

PrimaryKey::PrimaryKey(utility::vector1<Column> columns){
    columns_ = columns;
}

void PrimaryKey::add_column(Column column){
    columns_.push_back(column);
}

utility::vector1<Column> PrimaryKey::columns(){
    return this->columns_;
}
    
std::string PrimaryKey::print(){
    std::string primary_key_string = "PRIMARY KEY (";
    
    for(utility::vector1<Column>::const_iterator it=columns_.begin(); it != columns_.end(); ++it){
        if(it!=columns_.begin()){
            primary_key_string += ",";
        }
        primary_key_string += it->name();
    }
    
    primary_key_string += ")";
    return primary_key_string;
}
    
} // namespace sql_database
} // namespace utility


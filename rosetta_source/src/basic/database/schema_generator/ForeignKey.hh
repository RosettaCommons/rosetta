// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ForeignKey.hh
///
/// @brief
/// @author tim


#ifndef INCLUDED_basic_database_schema_generator_ForeignKey_HH
#define INCLUDED_basic_database_schema_generator_ForeignKey_HH

#include <basic/database/schema_generator/Column.hh>

//C++ Headers
#include <string>

namespace basic{
namespace database{
namespace schema_generator{
    
class ForeignKey
{
public:
    
    ForeignKey(Column column, std::string reference_table, std::string reference_column);
    
    ForeignKey(Column column, std::string reference_table, std::string reference_column, bool defer);
    
    void init_db_mode();
    
    std::string print();
    
    Column column();
    
private:
    
    std::string database_mode_;
    Column column_;
    std::string reference_column_;
    std::string reference_table_;
    bool defer_;
};

    
} // schema_generator
} // namespace database
} // namespace utility

#endif

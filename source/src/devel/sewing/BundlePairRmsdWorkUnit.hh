// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BundlePairRmsdWorkUnit.hh
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

#ifndef INCLUDED_devel_sewing_BundlePairRmsdWorkUnit_HH
#define INCLUDED_devel_sewing_BundlePairRmsdWorkUnit_HH

//Unit
#include <protocols/wum/DatabaseEntryWorkUnit.hh>
#include <protocols/wum/WorkUnitBase.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

//C++
#include <string>
#include <map>

namespace devel {
namespace sewing {

class BundlePairRmsdWorkUnit : public protocols::wum::DatabaseEntryWorkUnit {

public:
    
    BundlePairRmsdWorkUnit(utility::sql_database::sessionOP db_session);
    
    BundlePairRmsdWorkUnit( std::map<std::string,std::string> row_map ); 
    
    virtual ~BundlePairRmsdWorkUnit(){}
    
    virtual protocols::wum::WorkUnitBaseOP clone() const {
        return new BundlePairRmsdWorkUnit( *this );
    }
    
    /// @brief Calculate the pair RMSD using data from the results_map_
    virtual void run();
    
};

} //sewing namespace
} //devel namespace

#endif

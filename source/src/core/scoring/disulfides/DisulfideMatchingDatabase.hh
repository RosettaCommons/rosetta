// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/DisulfideMatchingDatabase.hh
/// @brief  Centroid Disulfide Energy Databases
/// @author rvernon@u.washington.edu
/// @date   02/09/10

#ifndef INCLUDED_core_scoring_disulfides_DisulfideMatchingDatabase_hh
#define INCLUDED_core_scoring_disulfides_DisulfideMatchingDatabase_hh

//Unit headers
//#include <core/scoring/disulfides/DisulfideMatchingDatabase.fwd.hh>

//Project headers
#include <core/types.hh>

//Utility Headers

#include <core/kinematics/RT.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace disulfides {

class DisulfideMatchingDatabase {
public:
	DisulfideMatchingDatabase();
	~DisulfideMatchingDatabase();

void
read_disulfide_database() const;
	//shouldn't really be const... move to setup for scoring later

utility::vector1< core::kinematics::RT > &
get_all_disulfides() const;

private:

	mutable bool db_init_;
	mutable utility::vector1< core::kinematics::RT > db_disulfides_;

}; // DisulfideMatchingDatabase

} // disulfides
} // scoring
} // core

#endif

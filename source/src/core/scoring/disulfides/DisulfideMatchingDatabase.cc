// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/CentroidDisulfideDatabase.cc
/// @brief  Centroid Disulfide Energy Databases
/// @author rvernon@u.washington.edu
/// @date   02/09/10

// Unit Headers
#include <core/scoring/disulfides/DisulfideMatchingDatabase.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <core/conformation/Atom.hh>
#include <core/kinematics/RT.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>


#include <basic/Tracer.hh>

#include <core/kinematics/Jump.hh>

static thread_local basic::Tracer TR( "core.scoring.disulfides.CentroidMatchingDatabase" );


using namespace core;
using core::scoring::disulfides::DisulfideMatchingDatabase;
using core::conformation::Residue;
using std::string;
using utility::vector1;

namespace core {
namespace scoring {
namespace disulfides {

/**
* Constructor
*/
DisulfideMatchingDatabase::DisulfideMatchingDatabase()
{
	db_init_ = false;
}

/**
* Deconstructor
*/
DisulfideMatchingDatabase::~DisulfideMatchingDatabase() {}

void
DisulfideMatchingDatabase::read_disulfide_database() const{

	utility::io::izstream disulfide_database;

	basic::database::open(disulfide_database, "sampling/disulfide_jump_database_wip.dat");

	std::string line;

	utility::vector1< core::kinematics::RT > disulfides;

	getline(disulfide_database, line);
	while ( !disulfide_database.eof() ) {
		std::istringstream line_stream(line);
		std::string junk;
		core::kinematics::RT rt;

		line_stream >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> rt;

		disulfides.push_back(rt);

		//std::cout << "RT " << rt << std::endl;

		getline(disulfide_database, line);
	}

	//std::cout << "READ DIS DB " << disulfides.size() << std::endl;

	db_disulfides_ = disulfides;
}

utility::vector1< core::kinematics::RT > &
DisulfideMatchingDatabase::get_all_disulfides() const {

	if ( db_init_ == false ) {
		//std::cout << "READING DB " << std::endl;
		read_disulfide_database();

		if ( db_disulfides_.size() == 0 ) {
			utility_exit_with_message("Failure to Read Disulfide Database");
		}

		db_init_ = true;
	}

	return db_disulfides_;
}

} // disulfides
} // scoring
} // core

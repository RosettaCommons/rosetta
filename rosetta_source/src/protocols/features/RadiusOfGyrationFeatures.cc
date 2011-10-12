// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/RadiusOfGyrationFeatures.cc
/// @brief  report the radius of gyration to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/RadiusOfGyrationFeatures.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/database/sql_utils.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols{
namespace features{

using std::string;
using core::Size;
using core::pose::Pose;
using core::scoring::methods::RG_Energy_Fast;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;

RadiusOfGyrationFeatures::RadiusOfGyrationFeatures(){}

RadiusOfGyrationFeatures::RadiusOfGyrationFeatures( RadiusOfGyrationFeatures const & ) :
	FeaturesReporter()
{}

RadiusOfGyrationFeatures::~RadiusOfGyrationFeatures(){}


string
RadiusOfGyrationFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS radius_of_gyration (\n"
		"	struct_id INTEGER,\n"
		"	radius_of_gyration REAL,\n"
		"	FOREIGN KEY(struct_id)\n"
		"		REFERENCES structures(struct_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id));";
}

Size
RadiusOfGyrationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	RG_Energy_Fast rg;

	statement stmt = (*db_session)
		<< "INSERT INTO radius_of_gyration VALUES (?,?);"
		<< struct_id
		<< rg.calculate_rg_score(pose, relevant_residues);
	basic::database::safely_write_to_database(stmt);
	return 0;
}

} // features namesapce
} // protocols namespace

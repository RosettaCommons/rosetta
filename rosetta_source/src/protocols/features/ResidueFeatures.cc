// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueFeatures.cc
/// @brief  report residue features to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueFeatures.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/database/sql_utils.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
//#include <cmath>

namespace protocols{
namespace features{

using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;


ResidueFeatures::ResidueFeatures() {}

ResidueFeatures::ResidueFeatures( ResidueFeatures const & src) :
	FeaturesReporter()
{}

ResidueFeatures::~ResidueFeatures()
{}

string
ResidueFeatures::type_name() const { return "ResidueFeatures"; }

string
ResidueFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS residues (\n"
		"	struct_id INTEGER,\n"
		"	resNum INTEGER,\n"
		"	name3 TEXT,\n"
		"	res_type TEXT,\n"
		"	FOREIGN KEY (struct_id)\n"
		"		REFERENCES structures (struct_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	CONSTRAINT resNum_is_positive CHECK (resNum >= 1),\n"
		"	PRIMARY KEY(struct_id, resNum));";
}

utility::vector1<std::string>
ResidueFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}

Size
ResidueFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	insert_residue_rows(pose, relevant_residues, struct_id, db_session);
	return 0;
}


void
ResidueFeatures::insert_residue_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){

	std::string statement_string = "INSERT INTO residues VALUES (?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
		Residue res = pose.residue(resNum);

		string const name3( res.name3() );
		string const res_type( res.name() );

		stmt.bind(1,struct_id);
		stmt.bind(2,resNum);
		stmt.bind(3,name3);
		stmt.bind(4,res_type);
		basic::database::safely_write_to_database(stmt);
	}
}

void
ResidueFeatures::delete_record(
	Size struct_id,
	sessionOP db_session) {

	delete_records_from_table("residues", struct_id, db_session);
}



} // namesapce
} // namespace

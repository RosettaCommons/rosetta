// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/GeometricSolvationFeatures.cc
/// @brief  report geometric solvation energy for each polar site to a features database
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/GeometricSolvationFeatures.hh>

//External

// Platform Headers
#include <core/pose/Pose.hh>
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>


// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>

namespace protocols {
namespace features {

using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::scoring::geometric_solvation::ExactOccludedHbondSolEnergy;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

GeometricSolvationFeatures::GeometricSolvationFeatures() :
	geo_sol_energy_(ExactOccludedHbondSolEnergy())
{

	//I would like to simply assert that this has been called, but that
	//currently is not possible.
	//geo_sol_energy_.setup_for_scoring(pose, scfxn_);

}

GeometricSolvationFeatures::GeometricSolvationFeatures(
	GeometricSolvationFeatures const & src ) :
	FeaturesReporter(),
	geo_sol_energy_(src.geo_sol_energy_)
{}

string
GeometricSolvationFeatures::type_name() const { return "GeometricSolvationFeatures"; }

void
GeometricSolvationFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_geometric_solvation_table_schema(db_session);
}

void
GeometricSolvationFeatures::write_geometric_solvation_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column hbond_site_id("hbond_site_id", DbDataTypeOP( new DbInteger() ));
	Column geometric_solvation_exact("geometric_solvation_exact", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(hbond_site_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(hbond_site_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("site_id");
	ForeignKey foreign_key(foreign_key_columns, "hbond_sites", reference_columns, true);

	Schema table("geometric_solvation", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(geometric_solvation_exact);

	table.write(db_session);
}

utility::vector1<std::string>
GeometricSolvationFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("HBondFeatures");
	return dependencies;
}

Size
GeometricSolvationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const &,
	StructureID struct_id,
	sessionOP db_session
){

	std::string select_string =
		"SELECT\n"
		"\tsite.site_id,\n"
		"\tsite.resNum,\n"
		"\tsite.atmNum\n"
		"FROM\n"
		"\thbond_sites as site\n"
		"WHERE\n"
		"\tsite.struct_id = ?;";
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	while ( res.next() ) {
		Size site_id, resNum, atmNum;
		res >> site_id >> resNum >> atmNum;

		Real const geometric_solvation_exact(
			geo_sol_energy_.compute_polar_group_sol_energy(
			pose,
			pose.residue(resNum),
			atmNum));

		std::string insert_string = "INSERT INTO geometric_solvation (struct_id, hbond_site_id, geometric_solvation_exact) VALUES (?,?,?)";
		statement insert_statement(basic::database::safely_prepare_statement(insert_string,db_session));
		insert_statement.bind(1,struct_id);
		insert_statement.bind(2,site_id);
		insert_statement.bind(3,geometric_solvation_exact);
		basic::database::safely_write_to_database(insert_statement);
	}

	// locate the polar sites from the hbond_sites table

	return 0;
}

} //namesapce
} //namespace

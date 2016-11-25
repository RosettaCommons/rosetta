// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProteinBackboneTorsionAngleFeatures.cc
/// @brief  report Backbone Torsional Angle features
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/ProteinBackboneTorsionAngleFeatures.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>


// Platform Headers
#include <core/pose/Pose.hh>

// External Headers
#include <cppdb/frontend.h>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/ProteinBackboneTorsionAngleFeaturesCreator.hh>

namespace protocols {
namespace features {

using std::string;
using cppdb::statement;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;


ProteinBackboneTorsionAngleFeatures::ProteinBackboneTorsionAngleFeatures(){}

ProteinBackboneTorsionAngleFeatures::ProteinBackboneTorsionAngleFeatures( ProteinBackboneTorsionAngleFeatures const & ) :
	FeaturesReporter()
{}

ProteinBackboneTorsionAngleFeatures::~ProteinBackboneTorsionAngleFeatures()= default;

// XRW TEMP string
// XRW TEMP ProteinBackboneTorsionAngleFeatures::type_name() const { return "ProteinBackboneTorsionAngleFeatures"; }

void
ProteinBackboneTorsionAngleFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_protein_backbone_torsion_angles_table_schema(db_session);
}

void
ProteinBackboneTorsionAngleFeatures::write_protein_backbone_torsion_angles_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum("resNum", DbDataTypeOP( new DbInteger() ));
	Column phi("phi", DbDataTypeOP( new DbReal() ));
	Column psi("psi", DbDataTypeOP( new DbReal() ));
	Column omega("omega", DbDataTypeOP( new DbReal() ));


	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(resNum);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("protein_backbone_torsion_angles", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(phi);
	table.add_column(psi);
	table.add_column(omega);

	table.write(db_session);
}

utility::vector1<std::string>
ProteinBackboneTorsionAngleFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
ProteinBackboneTorsionAngleFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	std::string statement_string ="INSERT INTO protein_backbone_torsion_angles (struct_id, resNum, phi, psi, omega) VALUES (?,?,?,?,?)";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( !check_relevant_residues( relevant_residues, i ) ) continue;

		Residue const & resi = pose.residue(i);
		if ( !resi.is_protein() ) continue;


		Real phi  (resi.mainchain_torsion(1));
		Real psi  (resi.mainchain_torsion(2));
		Real omega(resi.mainchain_torsion(3));

		stmt.bind(1,struct_id);
		stmt.bind(2,i);
		stmt.bind(3,phi);
		stmt.bind(4,psi);
		stmt.bind(5,omega);
		basic::database::safely_write_to_database(stmt);

	}
	return 0;
}

std::string ProteinBackboneTorsionAngleFeatures::type_name() const {
	return class_name();
}

std::string ProteinBackboneTorsionAngleFeatures::class_name() {
	return "ProteinBackboneTorsionAngleFeatures";
}

void ProteinBackboneTorsionAngleFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Reports the torsion angles (phi, psi, and omega) for a protein backbone, for every residue", attlist );
}

std::string ProteinBackboneTorsionAngleFeaturesCreator::type_name() const {
	return ProteinBackboneTorsionAngleFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ProteinBackboneTorsionAngleFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new ProteinBackboneTorsionAngleFeatures );
}

void ProteinBackboneTorsionAngleFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ProteinBackboneTorsionAngleFeatures::provide_xml_schema( xsd );
}


} // namesapce
} // namespace

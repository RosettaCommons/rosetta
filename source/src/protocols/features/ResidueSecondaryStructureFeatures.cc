// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ResidueSecondaryStructureFeatures.cc
/// @brief  report ResidueSecondaryStructure geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueSecondaryStructureFeatures.hh>

//External
#include <boost/assign/list_of.hpp>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/schema_generator/Constraint.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/ResidueSecondaryStructureFeaturesCreator.hh>

namespace protocols {
namespace features {

using std::string;
using core::scoring::dssp::Dssp;
using core::Size;
using core::conformation::Residue;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using basic::Tracer;

static Tracer TR("protocols.features.ResidueSecondaryStructureFeatures");

ResidueSecondaryStructureFeatures::ResidueSecondaryStructureFeatures() = default;

ResidueSecondaryStructureFeatures::ResidueSecondaryStructureFeatures(ResidueSecondaryStructureFeatures const &) :
	FeaturesReporter()
{}

ResidueSecondaryStructureFeatures::~ResidueSecondaryStructureFeatures() = default;

// XRW TEMP string
// XRW TEMP ResidueSecondaryStructureFeatures::type_name() const { return "ResidueSecondaryStructureFeatures"; }

void
ResidueSecondaryStructureFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;
	using namespace basic::database;
	using namespace boost::assign;

	Column code("code", DbDataTypeOP( new DbText(1) ), false);
	Column label("label", DbDataTypeOP( new DbText() ), false);
	Schema dssp_codes("dssp_codes", PrimaryKey(code));
	dssp_codes.add_column(label);

	dssp_codes.write(db_session);

	std::vector<std::string> dssp_cols = list_of("code")("label");
	insert_or_ignore("dssp_codes", dssp_cols, list_of("H")("H: a-Helix"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("E")("E: b-Sheet"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("T")("T: HB Turn"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("G")("G: 3/10 Helix"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("B")("B: b-Bridge"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("S")("S: Bend"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of("I")("I: pi-Helix"), db_session);
	insert_or_ignore("dssp_codes", dssp_cols, list_of(" ")("Irregular"), db_session);

	/******residue_secondary_structure******/
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column resNum("resNum", DbDataTypeOP( new DbInteger() ), false);
	Column dssp("dssp", DbDataTypeOP( new DbText(1) ));

	utility::vector1<Column> sec_struct_pkey_cols;
	sec_struct_pkey_cols.push_back(struct_id);
	sec_struct_pkey_cols.push_back(resNum);

	utility::vector1<Column> fkey_cols;
	fkey_cols.push_back(struct_id);
	fkey_cols.push_back(resNum);

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	ForeignKey dssp_fk(dssp, "dssp_codes", "code", true /*defer*/);

	Schema residue_secondary_structure("residue_secondary_structure", PrimaryKey(sec_struct_pkey_cols));

	residue_secondary_structure.add_column(struct_id);
	residue_secondary_structure.add_column(resNum);
	residue_secondary_structure.add_column(dssp);

	residue_secondary_structure.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));
	residue_secondary_structure.add_foreign_key(dssp_fk);

	residue_secondary_structure.write(db_session);
}

utility::vector1<std::string>
ResidueSecondaryStructureFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
ResidueSecondaryStructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	// compute dssp
	core::scoring::dssp::Dssp all_dssp(pose);

	//Create the statement strings outside the loops so we don't need to rcreate them for every residue
	std::string sec_structure_statement_string = "INSERT INTO residue_secondary_structure (struct_id, resNum, dssp) VALUES (?,?,?);";

	core::Size adjusted_resnum=0;
	for ( Size resNum=1; resNum <= pose.size(); ++resNum ) {
		if ( !check_relevant_residues( relevant_residues, resNum ) ) continue;


		//If this is not a protein residue then skip it. Keep a counter
		//of protein-only residues to reference DSSP
		if ( !pose.residue(resNum).is_protein() ) {
			continue;
		}
		adjusted_resnum++;

		string residue_secondary = string(1, all_dssp.get_dssp_secstruct(adjusted_resnum));
		statement stmt(basic::database::safely_prepare_statement(sec_structure_statement_string,db_session));
		stmt.bind(1,struct_id);
		stmt.bind(2,resNum);
		stmt.bind(3,residue_secondary);
		basic::database::safely_write_to_database(stmt);

	}
	return 0;
}

std::string ResidueSecondaryStructureFeatures::type_name() const {
	return class_name();
}

std::string ResidueSecondaryStructureFeatures::class_name() {
	return "ResidueSecondaryStructureFeatures";
}

void ResidueSecondaryStructureFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Report ResidueSecondaryStructure geometry and scores to features Statistics Scientific Benchmark",
		attlist );
}

std::string ResidueSecondaryStructureFeaturesCreator::type_name() const {
	return ResidueSecondaryStructureFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ResidueSecondaryStructureFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new ResidueSecondaryStructureFeatures );
}

void ResidueSecondaryStructureFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueSecondaryStructureFeatures::provide_xml_schema( xsd );
}



} // namesapce
} // namespace

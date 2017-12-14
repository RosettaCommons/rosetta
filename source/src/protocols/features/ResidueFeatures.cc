// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ResidueFeatures.cc
/// @brief  report residue features to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueFeatures.hh>

//External

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tools/make_vector.hh>

// External Headers
#include <cppdb/frontend.h>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/ResidueFeaturesCreator.hh>

// C++ Headers
//#include <cmath>

namespace protocols {
namespace features {

static basic::Tracer TR( "protocols.features.ResidueFeatures" );

using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

ResidueFeatures::ResidueFeatures() = default;

ResidueFeatures::ResidueFeatures( ResidueFeatures const & ) :
	FeaturesReporter()
{}

ResidueFeatures::~ResidueFeatures() = default;

// XRW TEMP string
// XRW TEMP ResidueFeatures::type_name() const { return "ResidueFeatures"; }

void
ResidueFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column resNum("resNum", DbDataTypeOP( new DbInteger() ), false);
	Column name3("name3", DbDataTypeOP( new DbText() ), false);
	Column res_type("res_type", DbDataTypeOP( new DbText() ), false);

	utility::vector1<Column> residues_pkey_cols;
	residues_pkey_cols.push_back(struct_id);
	residues_pkey_cols.push_back(resNum);

	Schema residues("residues", PrimaryKey(residues_pkey_cols));
	residues.add_column(struct_id);
	residues.add_column(resNum);
	residues.add_column(name3);
	residues.add_column(res_type);
	residues.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true));

	//TODO add constraint resNum > 0

	residues.write(db_session);

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
	StructureID const struct_id,
	sessionOP db_session
){
	insert_residue_rows(pose, relevant_residues, struct_id, db_session);
	return 0;
}


void
ResidueFeatures::insert_residue_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){


	InsertGenerator residues_insert("residues");
	residues_insert.add_column("struct_id");
	residues_insert.add_column("resNum");
	residues_insert.add_column("name3");
	residues_insert.add_column("res_type");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );

	for ( Size resNum=1; resNum <= pose.size(); ++resNum ) {
		if ( !check_relevant_residues(relevant_residues, resNum) ) continue;
		Residue res = pose.residue(resNum);

		string const name3( res.name3() );
		string const res_type( res.name() );

		RowDataBaseOP resnum_data( new RowData<Size>("resNum",resNum) );
		RowDataBaseOP name3_data( new RowData<string>("name3",name3) );
		RowDataBaseOP res_type_data( new RowData<string>("res_type",res_type) );

		residues_insert.add_row(utility::tools::make_vector(struct_id_data,resnum_data,name3_data,res_type_data));

	}

	residues_insert.write_to_database(db_session);
}

void
ResidueFeatures::delete_record(
	StructureID struct_id,
	sessionOP db_session) {

	delete_records_from_table("residues", struct_id, db_session);
}

std::string ResidueFeatures::type_name() const {
	return class_name();
}

std::string ResidueFeatures::class_name() {
	return "ResidueFeatures";
}

void ResidueFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Just record the residue number + type + name3 for each residue of a Pose", attlist );
}

std::string ResidueFeaturesCreator::type_name() const {
	return ResidueFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ResidueFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new ResidueFeatures );
}

void ResidueFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueFeatures::provide_xml_schema( xsd );
}



} // namesapce
} // namespace

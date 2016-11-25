// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/TrajectoryMapFeatures.cc
/// @brief  Map trajectory structure_ids to cycle number
/// @details Will only work properly when specially created by TrajectoryReportToDB
/// @author Kyle Barlow (kb@kylebarlow.com)

// Unit Headers
#include <protocols/features/TrajectoryMapFeatures.hh>

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
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/DbDataType.hh>
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
#include <protocols/features/TrajectoryMapFeaturesCreator.hh>

// C++ Headers
//#include <cmath>

namespace protocols {
namespace features {

static THREAD_LOCAL basic::Tracer TR( "protocols.features.TrajectoryMapFeatures" );

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

TrajectoryMapFeatures::TrajectoryMapFeatures() :
	current_cycle_(0)
{}

TrajectoryMapFeatures::TrajectoryMapFeatures( TrajectoryMapFeatures const & src) :
	FeaturesReporter(),
	current_cycle_(src.current_cycle_)
{}

TrajectoryMapFeatures::~TrajectoryMapFeatures() = default;

// XRW TEMP string
// XRW TEMP TrajectoryMapFeatures::type_name() const { return "TrajectoryMapFeatures"; }

void
TrajectoryMapFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column step("step", DbDataTypeOP( new DbInteger() ), false); // referred to as "cycle" in code

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(struct_id);
	pkey_cols.push_back(step);

	Schema trajectory_structures_steps("trajectory_structures_steps", PrimaryKey(pkey_cols));
	trajectory_structures_steps.add_column(struct_id);
	trajectory_structures_steps.add_column(step);
	trajectory_structures_steps.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true));

	trajectory_structures_steps.write(db_session);
}

utility::vector1<std::string>
TrajectoryMapFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}

Size
TrajectoryMapFeatures::report_features(
	Pose const &, // pose,
	vector1< bool > const &, // relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	InsertGenerator trajectory_structures_steps_insert("trajectory_structures_steps");
	trajectory_structures_steps_insert.add_column("struct_id");
	trajectory_structures_steps_insert.add_column("step");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id", struct_id) );
	RowDataBaseOP step_data( new RowData<Size>("step", current_cycle_) );

	trajectory_structures_steps_insert.add_row(
		utility::tools::make_vector(struct_id_data, step_data)
	);

	trajectory_structures_steps_insert.write_to_database(db_session);

	return 0;
}

void
TrajectoryMapFeatures::delete_record(
	StructureID struct_id,
	sessionOP db_session)
{
	delete_records_from_table("trajectory_structures_steps", struct_id, db_session);
}

void
TrajectoryMapFeatures::set_current_cycle(
	core::Size current_cycle
) {
	current_cycle_ = current_cycle;
}

core::Size
TrajectoryMapFeatures::get_current_cycle () const
{
	return current_cycle_;
}

std::string TrajectoryMapFeatures::type_name() const {
	return class_name();
}

std::string TrajectoryMapFeatures::class_name() {
	return "TrajectoryMapFeatures";
}

void TrajectoryMapFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Map trajectory structure_ids to cycle number", attlist );
}

std::string TrajectoryMapFeaturesCreator::type_name() const {
	return TrajectoryMapFeatures::class_name();
}

protocols::features::FeaturesReporterOP
TrajectoryMapFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new TrajectoryMapFeatures );
}

void TrajectoryMapFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TrajectoryMapFeatures::provide_xml_schema( xsd );
}


} // namesapce
} // namespace

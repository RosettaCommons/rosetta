// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/BetaTurnDetectionFeatures.cc
/// @brief report beta turns to a DB
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit Headers
#include <protocols/features/BetaTurnDetectionFeatures.hh>

// Package Headers
#include <protocols/features/BetaTurnDetection.hh>

// Project Headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/vector1.hh>

// Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// External Headers
#include <cppdb/frontend.h>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/BetaTurnDetectionFeaturesCreator.hh>

namespace protocols {
namespace features {

using std::string;
using basic::database::safely_write_to_database;
using basic::database::safely_prepare_statement;
using core::Size;
using core::SSize;
using core::pose::Pose;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;

BetaTurnDetectionFeatures::BetaTurnDetectionFeatures() : FeaturesReporter(), btd_( BetaTurnDetectionCOP( BetaTurnDetectionOP( new BetaTurnDetection ) ) ) {}

BetaTurnDetectionFeatures::BetaTurnDetectionFeatures( BetaTurnDetectionFeatures const & from ) :
	FeaturesReporter(), btd_( BetaTurnDetectionCOP( BetaTurnDetectionOP( new BetaTurnDetection( * from.btd_ ) ) ) ) {}

BetaTurnDetectionFeatures::~BetaTurnDetectionFeatures() = default;

// XRW TEMP string
// XRW TEMP BetaTurnDetectionFeatures::type_name() const { return "BetaTurnDetectionFeatures"; }

void
BetaTurnDetectionFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_beta_turns_table_schema( db_session );
}

void
BetaTurnDetectionFeatures::write_beta_turns_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id( "struct_id", DbDataTypeOP( new DbBigInt() ) );
	Column residue_begin( "residue_begin", DbDataTypeOP( new DbInteger() ) );
	Column turn_type( "turn_type", DbDataTypeOP( new DbText() ) );

	Columns primary_key_columns;
	primary_key_columns.push_back( struct_id );
	primary_key_columns.push_back( residue_begin );
	PrimaryKey primary_key( primary_key_columns );

	Columns foreign_key_columns;
	foreign_key_columns.push_back( struct_id );
	foreign_key_columns.push_back( residue_begin );
	vector1< string > reference_columns;
	reference_columns.push_back( "struct_id" );
	reference_columns.push_back( "resNum" );
	ForeignKey foreign_key( foreign_key_columns, "residues", reference_columns, true );

	Schema table( "beta_turns", primary_key );
	table.add_foreign_key( foreign_key );
	table.add_column( turn_type );

	table.write( db_session );
}

vector1< string >
BetaTurnDetectionFeatures::features_reporter_dependencies() const {
	vector1< string > dependencies;
	dependencies.push_back( "ResidueFeatures" );
	return dependencies;
}

/// @details
/// An anchor is a take off and landing for a loop.
/// Every residue in the loop must be relevant in order for the loop to be stored.
Size
BetaTurnDetectionFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	string beta_turns_stmt_string = "INSERT INTO beta_turns (struct_id, residue_begin, turn_type) VALUES (?,?,?);";
	statement beta_turns_stmt(
		safely_prepare_statement( beta_turns_stmt_string, db_session ) );

	for ( SSize begin = 1; begin <= SSize( pose.size() - btd_->beta_turn_length() ); ++begin ) {
		Size end = begin + btd_->beta_turn_length();

		if ( ! check_relevant_residues_range( relevant_residues, begin, end ) ||
				! btd_->residue_range_is_protein( pose, begin, end ) ||
				! btd_->all_turn_residues_are_on_the_same_chain( pose, begin ) ||
				! btd_->beta_turn_present( pose, begin ) ) {
			continue;
		}

		// Add stuff to database
		beta_turns_stmt.bind( 1, struct_id );
		beta_turns_stmt.bind( 2, begin );
		beta_turns_stmt.bind( 3, btd_->beta_turn_type( pose, begin ) );
		safely_write_to_database( beta_turns_stmt );

	}
	return 0;
}

std::string BetaTurnDetectionFeatures::type_name() const {
	return class_name();
}

std::string BetaTurnDetectionFeatures::class_name() {
	return "BetaTurnDetectionFeatures";
}

void BetaTurnDetectionFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Gives features about beta turn geometry", attlist );
}

std::string BetaTurnDetectionFeaturesCreator::type_name() const {
	return BetaTurnDetectionFeatures::class_name();
}

protocols::features::FeaturesReporterOP
BetaTurnDetectionFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new BetaTurnDetectionFeatures );
}

void BetaTurnDetectionFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BetaTurnDetectionFeatures::provide_xml_schema( xsd );
}


} // namesapce features
} // namespace protocols

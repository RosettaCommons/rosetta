// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>
#include <protocols/features/FeaturesReporter.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

--class--::--class--():
 	protocols::features::FeaturesReporter()
{

}

--class--::~--class--() {}

--class--::--class--( --class-- const & src ) :
	protocols::features::FeaturesReporter( src )
{}

--class--OP
--class--::clone() const {
	return --class--OP( new --class--( *this ) );
}

void
--class--::write_schema_to_db( utility::sql_database::sessionOP db_session ) const {
	// implement me!
}

utility::vector1< std::string >
--class--::features_reporter_dependencies() const {

}

Size
--class--::report_features(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & relevant_residues,
	protocols::features::StructureID const struct_id,
	utility::sql_database::sessionOP db_session
) {

}

std::string --class--::type_name() const {
	return class_name();
}

std::string --class--::class_name() {
	return "--class--";
}

void --class--::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}

// Methods for the features reporter creator

--class--Creator::--class--Creator() {}

--class--Creator::~--class--Creator() {}

FeaturesReporterOP
--class--Creator::create_features_reporter() const {
	return FeaturesReporterOP( new --class-- );
}

std::string
--class--Creator::type_name() const {
	return --class--::class_name();
}

void --class--Creator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	--class--::provide_xml_schema( xsd );
}


--end_namespace--

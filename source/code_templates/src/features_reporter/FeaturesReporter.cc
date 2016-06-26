// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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

static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );

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
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
) {
	
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
	return "--class--";
}

--end_namespace--




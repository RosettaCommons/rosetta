// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file features_database_schema.cc
/// @brief simple application for interconverting Rosetta recognized structure formats
/// @author Matthew O'Meara


// Project Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/sql_utils.hh>

// Platform Headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/FeaturesReporterFactory.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <core/types.hh>

using std::string;
using std::cout;
using std::endl;
using core::Size;
using utility::sql_database::sessionOP;
using utility::vector1;
using basic::database::get_db_session;
using protocols::features::FeaturesReporterOP;
using protocols::features::FeaturesReporterFactory;

void
write_schema_to_db(
	sessionOP db_session
) {
	FeaturesReporterFactory * features_reporter_factory(
		FeaturesReporterFactory::get_instance());

	vector1<string> features_reporter_names(
		features_reporter_factory->get_all_features_names());
	for(Size i=1; i <= features_reporter_names.size(); ++i){
		FeaturesReporterOP features_reporter(
			features_reporter_factory->get_features_reporter(
				features_reporter_names[i]));

		features_reporter->write_schema_to_db(db_session);
	}
}

int
main( int argc, char* argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// define relevant options

	try {
		// options, random initialization
		devel::init( argc, argv );


		sessionOP db_session(get_db_session());
		write_schema_to_db(db_session);

	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
    return 0;
} // main

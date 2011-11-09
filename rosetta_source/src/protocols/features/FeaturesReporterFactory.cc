// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/FeaturesReporterFactory.cc
/// @brief  Factory for creating FeaturesReporters objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/FeaturesReporterFactory.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/FeaturesReporterCreator.hh>

// Package Headers
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <sstream>






namespace protocols {
namespace features {

using std::endl;
using std::string;
using std::stringstream;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using protocols::filters::Filters_map;
using protocols::moves::DataMap;
using protocols::moves::Movers_map;
using utility::tag::TagPtr;

static basic::Tracer tr("protocols.features.FeaturesReporterFactory");

FeaturesReporterFactory * FeaturesReporterFactory::instance_( 0 );

/// @details Private constructor insures correctness of singleton.
FeaturesReporterFactory::FeaturesReporterFactory() {}

FeaturesReporterFactory::FeaturesReporterFactory(
	const FeaturesReporterFactory &
) {}

FeaturesReporterFactory::~FeaturesReporterFactory() {}


FeaturesReporterFactory *
FeaturesReporterFactory::get_instance()
{
	if ( instance_ == 0 ) {
		instance_ = new FeaturesReporterFactory;
	}
	return instance_;
}


void
FeaturesReporterFactory::factory_register(
	FeaturesReporterCreatorCOP creator
) {
	types_[ creator->type_name() ] = creator;
}


FeaturesReporterOP
FeaturesReporterFactory::get_features_reporter(
	string const & type_name
) {
	tr.Trace << "generate features reporter of type " << type_name << std::endl;
	FeaturesReporterCreatorMap::const_iterator iter = types_.find( type_name );
	if (iter != types_.end()) {
		return iter->second->create_features_reporter();
	} else {

		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized FeaturesReporter "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new FeaturesReporter in the FeaturesReporterFactory" << endl
			<< "known FeaturesReporter types are:" << endl;

		foreach(const FeaturesReporterCreatorMap::value_type& type, types_){
			error_msg << "\t" << type.first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}

FeaturesReporterOP
FeaturesReporterFactory::get_features_reporter(
	TagPtr const tag,
	DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose
) {
	assert(tag->getName() == "feature");

	string type_name;
	if(!tag->hasOption("name")){
		utility_exit_with_message("'feature' tags require a name field");
	} else {
		type_name = tag->getOption<string>("name");
	}

	FeaturesReporterOP features_reporter(get_features_reporter(type_name));

	features_reporter->parse_my_tag(tag, data, filters, movers, pose);
	return features_reporter;
}

} // namespace
} // namespace

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
// AUTO-REMOVED #include <protocols/moves/Mover.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ Headers
#include <sstream>

//Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH




namespace protocols {
namespace features {

using std::endl;
using std::string;
using std::stringstream;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using utility::tag::TagCOP;

static basic::Tracer tr("protocols.features.FeaturesReporterFactory");

FeaturesReporterFactory * FeaturesReporterFactory::instance_( 0 );

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex FeaturesReporterFactory::singleton_mutex_;

std::mutex & FeaturesReporterFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
FeaturesReporterFactory * FeaturesReporterFactory::get_instance()
{
	boost::function< FeaturesReporterFactory * () > creator = boost::bind( &FeaturesReporterFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

FeaturesReporterFactory *
FeaturesReporterFactory::create_singleton_instance()
{
	return new FeaturesReporterFactory;
}

/// @details Private constructor insures correctness of singleton.
FeaturesReporterFactory::FeaturesReporterFactory() {}

FeaturesReporterFactory::FeaturesReporterFactory(
	const FeaturesReporterFactory &
) {}

FeaturesReporterFactory::~FeaturesReporterFactory() {}

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
utility::vector1<std::string> FeaturesReporterFactory::get_all_features_names()
{
	utility::vector1<std::string> collection;
	FeaturesReporterCreatorMap::const_iterator iter = types_.begin();
	while ( iter != types_.end() ) {
		collection.push_back(iter->first);
		iter++;
	}
	return collection;

}
FeaturesReporterOP
FeaturesReporterFactory::get_features_reporter(
	TagCOP const tag,
	basic::datacache::DataMap & data,
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

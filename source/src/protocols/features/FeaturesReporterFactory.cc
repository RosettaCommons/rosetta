// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/FeaturesReporterFactory.cc
/// @brief  Factory for creating FeaturesReporters objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/FeaturesReporterFactory.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/FeaturesReporterCreator.hh>
#include <protocols/features/feature_schemas.hh>
// Package Headers
#include <basic/Tracer.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ Headers
#include <sstream>

//Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>


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

static basic::Tracer tr( "protocols.features.FeaturesReporterFactory" );

/// @details Private constructor insures correctness of singleton.
FeaturesReporterFactory::FeaturesReporterFactory() {}

FeaturesReporterFactory::~FeaturesReporterFactory() = default;

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
	if ( iter != types_.end() ) {
		return iter->second->create_features_reporter();
	} else {
		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized FeaturesReporter "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new FeaturesReporter in the FeaturesReporterFactory" << endl
			<< "known FeaturesReporter types are:" << endl;

		for ( auto const & type : types_ ) {
			error_msg << "\t" << type.first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return nullptr;
}
utility::vector1<std::string> FeaturesReporterFactory::get_all_features_names()
{
	utility::vector1<std::string> collection;
	FeaturesReporterCreatorMap::const_iterator iter = types_.begin(), end = types_.end();
	while ( iter != end ) {
		collection.push_back(iter->first);
		++iter;
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
	//debug_assert(tag->getName() == "feature");

	string type_name = tag->getName();
	/*
	if ( !tag->hasOption("name") ) {
	utility_exit_with_message("'feature' tags require a name field");
	} else {
	type_name = tag->getOption<string>("name");
	}
	*/
	FeaturesReporterOP features_reporter(get_features_reporter(type_name));

	features_reporter->parse_my_tag(tag, data, filters, movers, pose);
	return features_reporter;
}
void
FeaturesReporterFactory::define_features_reporter_xml_schema_group( utility::tag::XMLSchemaDefinition & xsd ) const{
	try{
		utility::tag::define_xml_schema_group(
			types_,
			features_reporter_xml_schema_group_name(),
			& complex_type_name_for_features_reporter,
			xsd );
	} catch( utility::excn::EXCN_Msg_Exception const & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Could not generate an XML Schema for FeaturesReporter from FeaturesReporterFactory; offending class"
			" must call protocols::features::complex_type_name_for_features_reporter when defining"
			" its XML Schema\n" + e.msg() );
	}
}

std::string
FeaturesReporterFactory::features_reporter_xml_schema_group_name(){
	return "features_reporter";
}


} // namespace
} // namespace

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/PosePropertyReporterFactory.cc
/// @brief Factory for PosePropertyReporters
/// @author Luki Goldschmidt <lugo@uw.edu>

// Unit Headers
#include <protocols/rosetta_scripts/PosePropertyReporterFactory.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.fwd.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR( "protocols.rosetta_scripts.PosePropertyReporterFactory" );

PosePropertyReporterFactory::PosePropertyReporterFactory(){}

PosePropertyReporterFactory::~PosePropertyReporterFactory()= default;

/// @brief add a PosePropertyReporter prototype, using its default type name as the map key
void
PosePropertyReporterFactory::factory_register( PosePropertyReporterCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	std::string const pose_selector_type( creator->keyname() );
	if ( pose_selector_type == "UNDEFINED NAME" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Can't map derived PosePropertyReporter with undefined type name.");
	}
	if ( reporter_creator_map_.find( pose_selector_type ) != reporter_creator_map_.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption("PosePropertyReporterFactory::factory_register already has a pose selector creator with name \"" + pose_selector_type + "\".  Conflicting pose selector names" );
	}
	reporter_creator_map_[ pose_selector_type ] = creator;
}


/// @brief return new PosePropertyReporter by key lookup in reporter_creator_map_ (new PosePropertyReporter parses Tag if provided)
PosePropertyReporterOP
PosePropertyReporterFactory::newPosePropertyReporter( std::string const & pose_selector_type )
{
	PosePropertyReporterMap::const_iterator iter( reporter_creator_map_.find( pose_selector_type ) );
	if ( iter != reporter_creator_map_.end() ) {
		if ( ! iter->second ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: PosePropertyReporterCreatorOP prototype for " + pose_selector_type + " is NULL!" );
		}
		return iter->second->create_reporter();
	} else {
		TR<<"Available pose property reporters: ";
		for ( PosePropertyReporterMap::const_iterator it = reporter_creator_map_.begin(); it != reporter_creator_map_.end(); ++it ) {
			TR<<it->first<<", ";
		}
		TR<<std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( pose_selector_type + " is not known to the PosePropertyReporterFactory. Was it registered via a PosePropertyReporterRegistrator in one of the init.cc files?" );
		return nullptr;
	}
}

/// @brief return new PosePropertyReporter by Tag parsing
PosePropertyReporterOP
PosePropertyReporterFactory::newPosePropertyReporter(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {
	PosePropertyReporterOP selector( newPosePropertyReporter( tag->getName() ) );
	runtime_assert( selector != nullptr );
	selector->parse_my_tag( tag, data, filters, movers, pose );
	return selector;
}



void
PosePropertyReporterFactory::define_pose_reporter_group( utility::tag::XMLSchemaDefinition & xsd ) const{
	try{
		utility::tag::define_xml_schema_group(
			reporter_creator_map_,
			pose_reporter_group_name(),
			& complex_type_name_for_pose_reporter,
			xsd );
	} catch( utility::excn::EXCN_Msg_Exception const & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Could not generate an XML Schema for PosePropertyReporter from PosePropertyReporterFactory; offending class"
			" must call protocols::rosetta_scripts::complex_type_name_for_pose_reporter when defining"
			" its XML Schema\n" + e.msg() );
	}
}

std::string
PosePropertyReporterFactory::pose_reporter_group_name(){
	return "pose_property_reporter";
}

std::string
PosePropertyReporterFactory::complex_type_name_for_pose_reporter( std::string const & reporter_name ){
	return "pose_property_reporter_" + reporter_name + "_complex_type";
}





} //namespace rosetta_scripts
} //namespace protocols

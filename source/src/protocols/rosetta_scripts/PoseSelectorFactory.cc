// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/PoseSelectorFactory.cc
/// @brief  Factory for PoseSelectors
/// @author Luki Goldschmidt <lugo@uw.edu>

// Unit Headers
#include <protocols/rosetta_scripts/PoseSelectorFactory.hh>
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>
#include <protocols/rosetta_scripts/PoseSelector.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR( "protocols.rosetta_scripts.PoseSelectorFactory" );

PoseSelectorFactory::PoseSelectorFactory(){}

PoseSelectorFactory::~PoseSelectorFactory()= default;

/// @brief add a PoseSelector prototype, using its default type name as the map key
void
PoseSelectorFactory::factory_register( PoseSelectorCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	std::string const pose_selector_type( creator->keyname() );
	if ( pose_selector_type == "UNDEFINED NAME" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Can't map derived PoseSelector with undefined type name.");
	}
	if ( poseselector_creator_map_.find( pose_selector_type ) != poseselector_creator_map_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "PoseSelectorFactory::factory_register already has a pose selector creator with name \"" + pose_selector_type + "\".  Conflicting pose selector names" );
	}
	poseselector_creator_map_[ pose_selector_type ] = creator;
}


/// @brief return new PoseSelector by key lookup in poseselector_creator_map_ (new PoseSelector parses Tag if provided)
PoseSelectorOP
PoseSelectorFactory::newPoseSelector( std::string const & pose_selector_type )
{
	PoseSelectorMap::const_iterator iter( poseselector_creator_map_.find( pose_selector_type ) );
	if ( iter != poseselector_creator_map_.end() ) {
		if ( ! iter->second ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Error: PoseSelectorCreatorOP prototype for " + pose_selector_type + " is NULL!" );
		}
		return iter->second->create_selector();
	} else {
		TR<<"Available pose selectors: ";
		for ( PoseSelectorMap::const_iterator it = poseselector_creator_map_.begin(); it != poseselector_creator_map_.end(); ++it ) {
			TR<<it->first<<", ";
		}
		TR<<std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  pose_selector_type + " is not known to the PoseSelectorFactory. Was it registered via a PoseSelectorRegistrator in one of the init.cc files?" );
		return nullptr;
	}
}

/// @brief return new PoseSelector by Tag parsing
PoseSelectorOP
PoseSelectorFactory::newPoseSelector(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {
	PoseSelectorOP selector( newPoseSelector( tag->getName() ) );
	runtime_assert( selector != nullptr );
	selector->parse_my_tag( tag, data, filters, movers, pose );
	return selector;
}



void
PoseSelectorFactory::define_pose_selector_group( utility::tag::XMLSchemaDefinition & xsd ) const{
	try{
		utility::tag::define_xml_schema_group(
			poseselector_creator_map_,
			pose_selector_group_name(),
			& complex_type_name_for_pose_selector,
			xsd );
	} catch( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for PoseSelector from PoseSelectorFactory; offending class"
			" must call protocols::rosetta_scripts::complex_type_name_for_pose_selector when defining"
			" its XML Schema\n" + e.msg() );
	}
}

std::string
PoseSelectorFactory::pose_selector_group_name(){
	return "pose_selector";
}
std::string
PoseSelectorFactory::complex_type_name_for_pose_selector( std::string const & selector_name ){
	return "pose_selector_" + selector_name + "_complex_type";
}


} //namespace rosetta_scripts
} //namespace protocols

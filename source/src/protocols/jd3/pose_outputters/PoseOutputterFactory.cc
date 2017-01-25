// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_outputters/PoseOutputterFactory.cc
/// @brief  PoseOutputterFactory class that holds the list of classes able to create Poses from inputs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>

#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputterCreator.hh>

// required for passing to PoseOutputter::parse_my_tag
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, throw utility::excn::EXCN_RosettaScriptsOption
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
#include <utility/excn/Exceptions.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>


namespace protocols {
namespace jd3 {
namespace pose_outputters {


static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.pose_outputters.PoseOutputterFactory" );

PoseOutputterFactory::PoseOutputterFactory() :
	throw_on_double_registration_( false )
{
}

//PoseOutputterFactory::~PoseOutputterFactory(){}

/// @brief add a PoseOutputter prototype, using its default type name as the map key
void
PoseOutputterFactory::factory_register( PoseOutputterCreatorOP creator )
{
	runtime_assert( creator != 0 );
	std::string const pose_outputter_type( creator->keyname() );
	if ( pose_outputter_creator_map_.find( pose_outputter_type ) != pose_outputter_creator_map_.end() ) {
		std::string err_msg = "PoseOutputterFactory::factory_register already has a pose_outputter creator with name \""
			+ pose_outputter_type + "\".  Conflicting PoseOutputter names";
		if ( throw_on_double_registration_ ) {
			throw utility::excn::EXCN_Msg_Exception( err_msg );
		} else {
			utility_exit_with_message( err_msg );
		}
	}
	pose_outputter_creator_map_[ pose_outputter_type ] = creator;
	creator_list_.push_back( creator );
}

void
PoseOutputterFactory::factory_register( SecondaryPoseOutputterCreatorOP creator )
{
	runtime_assert( creator != 0 );
	std::string const secondary_pose_outputter_type( creator->keyname() );
	if ( secondary_pose_outputter_creator_map_.find( secondary_pose_outputter_type ) != secondary_pose_outputter_creator_map_.end() ) {
		std::string err_msg = "PoseOutputterFactory::factory_register already has a secondary pose outputter creator with name \""
			+ secondary_pose_outputter_type + "\".  Conflicting SecondaryPoseOutputter names";
		if ( throw_on_double_registration_ ) {
			throw utility::excn::EXCN_Msg_Exception( err_msg );
		} else {
			utility_exit_with_message( err_msg );
		}
	}
	secondary_pose_outputter_creator_map_[ secondary_pose_outputter_type ] = creator;
	secondary_creator_list_.push_back( creator );
}


/// @details return new PoseOutputter by key lookup in pose_outputter_prototype_map_
PoseOutputterOP
PoseOutputterFactory::new_pose_outputter( std::string const & pose_outputter_type ) const
{
	PoseOutputterMap::const_iterator iter( pose_outputter_creator_map_.find( pose_outputter_type ) );
	if ( iter != pose_outputter_creator_map_.end() ) {
		if ( ! iter->second ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: PoseOutputterCreatorOP prototype for " + pose_outputter_type + " is NULL!" );
		}
		return iter->second->create_outputter();
	} else {
		TR<<"Available PoseOutputters: ";
		for ( PoseOutputterMap::const_iterator pose_outputter_it = pose_outputter_creator_map_.begin();
				pose_outputter_it != pose_outputter_creator_map_.end(); ++pose_outputter_it ) {
			TR << pose_outputter_it->first << ", ";
		}
		TR << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( pose_outputter_type + " is not known to the PoseOutputterFactory."
			" Was it registered via a PoseOutputterRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return PoseOutputterOP();
	}
}

/// @details return new SecondaryPoseOutputter by key lookup in secondary_pose_outputter_prototype_map_
SecondaryPoseOutputterOP
PoseOutputterFactory::new_secondary_outputter( std::string const & secondary_outputter_type ) const
{
	SecondaryOutputterMap::const_iterator iter( secondary_pose_outputter_creator_map_.find( secondary_outputter_type ) );
	if ( iter != secondary_pose_outputter_creator_map_.end() ) {
		if ( ! iter->second ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: SecondaryPoseOutputterCreatorOP prototype for " + secondary_outputter_type + " is NULL!" );
		}
		return iter->second->create_outputter();
	} else {
		TR<<"Available SecondaryPoseOutputters: ";
		for ( SecondaryOutputterMap::const_iterator secondary_pose_outputter_it = secondary_pose_outputter_creator_map_.begin();
				secondary_pose_outputter_it != secondary_pose_outputter_creator_map_.end(); ++secondary_pose_outputter_it ) {
			TR << secondary_pose_outputter_it->first << ", ";
		}
		TR << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( secondary_outputter_type + " is not known to the PoseOutputterFactory."
			" Was it registered via a SecondaryPoseOutputterRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return SecondaryPoseOutputterOP();
	}
}

/// @details In the event that there are multiple PoseOutputters specified on the command line,
/// the one that is returned is the first in alphabetical order by key that is specified.
PoseOutputterOP
PoseOutputterFactory::pose_outputter_from_command_line() const
{
	for ( PoseOutputterMap::const_iterator pose_outputter_it = pose_outputter_creator_map_.begin();
			pose_outputter_it != pose_outputter_creator_map_.end(); ++pose_outputter_it ) {
		if ( pose_outputter_it->second->outputter_specified_by_command_line() ) {
			return pose_outputter_it->second->create_outputter();
		}
	}
	// having reached here, no pose outputter is explicity set on the command line, and so
	// the default outputter -- the PDBPoseOutputter -- should be returned
	return PoseOutputterOP( new PDBPoseOutputter );
}

/// @details In the event that there are multiple PoseOutputters specified on the command line,
/// the one that is returned is the first in alphabetical order by key that is specified.
std::list< SecondaryPoseOutputterOP >
PoseOutputterFactory::secondary_pose_outputters_from_command_line() const
{
	std::list< SecondaryPoseOutputterOP > secondary_outputters;
	for ( SecondaryOutputterMap::const_iterator secondary_pose_outputter_it = secondary_pose_outputter_creator_map_.begin();
			secondary_pose_outputter_it != secondary_pose_outputter_creator_map_.end(); ++secondary_pose_outputter_it ) {
		if ( secondary_pose_outputter_it->second->outputter_specified_by_command_line() ) {
			secondary_outputters.push_back( secondary_pose_outputter_it->second->create_outputter() );
		}
	}
	return secondary_outputters;
}


void PoseOutputterFactory::set_throw_on_double_registration()
{
	throw_on_double_registration_ = true;
}


void PoseOutputterFactory::define_pose_outputter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	utility::tag::define_xml_schema_group(
		pose_outputter_creator_map_,
		pose_outputter_xml_schema_group_name(),
		& PoseOutputterFactory::complex_type_name_for_pose_outputter,
		xsd );

	utility::tag::define_xml_schema_group(
		secondary_pose_outputter_creator_map_,
		secondary_pose_outputter_xml_schema_group_name(),
		& PoseOutputterFactory::complex_type_name_for_secondary_pose_outputter,
		xsd );

}

void PoseOutputterFactory::define_secondary_pose_outputter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	utility::tag::define_xml_schema_group(
		secondary_pose_outputter_creator_map_,
		secondary_pose_outputter_xml_schema_group_name(),
		& PoseOutputterFactory::complex_type_name_for_secondary_pose_outputter,
		xsd );
}

void PoseOutputterFactory::list_outputter_options_read( utility::options::OptionKeyList & read_options ) const
{
	for ( CreatorList::const_iterator iter = creator_list_.begin(); iter != creator_list_.end(); ++iter ) {
		(*iter)->list_options_read( read_options );
	}
}

void PoseOutputterFactory::list_secondary_outputter_options_read( utility::options::OptionKeyList & read_options ) const
{
	for ( SecondaryCreatorList::const_iterator iter = secondary_creator_list_.begin();
			iter != secondary_creator_list_.end(); ++iter ) {
		(*iter)->list_options_read( read_options );
	}
}


PoseOutputterFactory::CreatorList const &
PoseOutputterFactory::pose_outputter_creators() const
{
	return creator_list_;
}

PoseOutputterFactory::SecondaryCreatorList const &
PoseOutputterFactory::secondary_pose_outputter_creators() const
{
	return secondary_creator_list_;
}


std::string PoseOutputterFactory::pose_outputter_xml_schema_group_name()
{
	return "pose_outputter";
}

std::string PoseOutputterFactory::complex_type_name_for_pose_outputter( std::string const & outputter_key )
{
	return "pose_outputter_" + outputter_key + "_type";
}

std::string PoseOutputterFactory::secondary_pose_outputter_xml_schema_group_name()
{
	return "secondary_pose_outputter";
}

std::string PoseOutputterFactory::complex_type_name_for_secondary_pose_outputter( std::string const & outputter_key )
{
	return "secondary_pose_outputter_" + outputter_key + "_type";
}


} //namespace pose_outputters
} //namespace jd3
} //namespace protocols

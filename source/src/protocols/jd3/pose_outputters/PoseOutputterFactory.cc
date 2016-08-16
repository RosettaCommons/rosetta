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

// required for passing to PoseOutputter::parse_my_tag
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, throw utility::excn::EXCN_RosettaScriptsOption
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>


namespace protocols {
namespace jd3 {
namespace pose_outputters {


static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.pose_outputters.PoseOutputterFactory" );

#if defined MULTI_THREADED && defined CXX11
std::atomic< PoseOutputterFactory * > PoseOutputterFactory::instance_( 0 );
#else
PoseOutputterFactory * PoseOutputterFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex PoseOutputterFactory::singleton_mutex_;

std::mutex & PoseOutputterFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of (pointer to) this singleton class
PoseOutputterFactory * PoseOutputterFactory::get_instance()
{
	boost::function< PoseOutputterFactory * () > creator = boost::bind( &PoseOutputterFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

PoseOutputterFactory *
PoseOutputterFactory::create_singleton_instance()
{
	return new PoseOutputterFactory;
}

PoseOutputterFactory::PoseOutputterFactory() :
	throw_on_double_registration_( false )
{
}

PoseOutputterFactory::~PoseOutputterFactory(){}

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


/// @brief return new PoseOutputter by key lookup in pose_outputter_prototype_map_ (new PoseOutputter parses Tag if provided)
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

/// @details In the event that there are multiple PoseOutputters specified on the command line,
/// the one that is returned is the first in alphabetical order by key that is specified.
PoseOutputterOP
PoseOutputterFactory::pose_outputter_from_command_line() const
{
	for ( PoseOutputterMap::const_iterator pose_outputter_it = pose_outputter_creator_map_.begin();
			pose_outputter_it != pose_outputter_creator_map_.end(); ++pose_outputter_it ) {
		PoseOutputterOP outputter = pose_outputter_it->second->create_outputter();
		if ( outputter->outputter_specified_by_command_line() ) {
			return outputter;
		}
	}
	// having reached here, no pose outputter is explicity set on the command line, and so
	// the default outputter -- the PDBPoseOutputter -- should be returned
	return PoseOutputterOP( new PDBPoseOutputter );
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
}

void PoseOutputterFactory::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	for ( CreatorList::const_iterator iter = creator_list_.begin(); iter != creator_list_.end(); ++iter ) {
		(*iter)->list_options_read( read_options );
	}
}


std::string PoseOutputterFactory::pose_outputter_xml_schema_group_name()
{
	return "pose_outputter";
}

std::string PoseOutputterFactory::complex_type_name_for_pose_outputter( std::string const & outputter_key )
{
	return "pose_outputter_" + outputter_key + "_type";
}

} //namespace pose_outputters
} //namespace jd3
} //namespace protocols

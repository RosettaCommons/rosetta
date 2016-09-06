// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/PoseInputterFactory.cc
/// @brief  PoseInputterFactory class that holds the list of classes able to create Poses from inputs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <protocols/jd3/pose_inputters/PoseInputterFactory.hh>

#include <protocols/jd3/pose_inputters/PoseInputter.hh>
#include <protocols/jd3/pose_inputters/PoseInputterCreator.hh>

// core headers
#include <core/types.hh>

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
namespace pose_inputters {


static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.pose_inputters.PoseInputterFactory" );

#ifdef MULTI_THREADED
std::atomic< PoseInputterFactory * > PoseInputterFactory::instance_( 0 );
#else
PoseInputterFactory * PoseInputterFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED

std::mutex PoseInputterFactory::singleton_mutex_;

std::mutex & PoseInputterFactory::singleton_mutex() { return singleton_mutex_; }

#endif

/// @brief static function to get the instance of (pointer to) this singleton class
PoseInputterFactory * PoseInputterFactory::get_instance()
{
	boost::function< PoseInputterFactory * () > creator = boost::bind( &PoseInputterFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

PoseInputterFactory *
PoseInputterFactory::create_singleton_instance()
{
	return new PoseInputterFactory;
}

PoseInputterFactory::PoseInputterFactory() :
	throw_on_double_registration_( false )
{
}

PoseInputterFactory::~PoseInputterFactory(){}

/// @brief add a PoseInputter prototype, using its default type name as the map key
void
PoseInputterFactory::factory_register( PoseInputterCreatorOP creator )
{
	runtime_assert( creator != 0 );
	std::string const pose_inputter_type( creator->keyname() );
	if ( pose_inputter_creator_map_.find( pose_inputter_type ) != pose_inputter_creator_map_.end() ) {
		std::string err_msg = "PoseInputterFactory::factory_register already has a pose_inputter creator with name \""
			+ pose_inputter_type + "\".  Conflicting PoseInputter names";
		if ( throw_on_double_registration_ ) {
			throw utility::excn::EXCN_Msg_Exception( err_msg );
		} else {
			utility_exit_with_message( err_msg );
		}
	}
	pose_inputter_creator_map_[ pose_inputter_type ] = creator;
	creator_list_.push_back( creator );
}


/// @brief return new PoseInputter by key lookup in pose_inputter_prototype_map_ (new PoseInputter parses Tag if provided)
PoseInputterOP
PoseInputterFactory::new_pose_inputter( std::string const & pose_inputter_type ) const
{
	PoseInputterMap::const_iterator iter( pose_inputter_creator_map_.find( pose_inputter_type ) );
	if ( iter != pose_inputter_creator_map_.end() ) {
		if ( ! iter->second ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: PoseInputterCreatorOP prototype for " + pose_inputter_type + " is NULL!" );
		}
		return iter->second->create_inputter();
	} else {
		TR<<"Available PoseInputters: ";
		for ( PoseInputterMap::const_iterator pose_inputter_it = pose_inputter_creator_map_.begin();
				pose_inputter_it != pose_inputter_creator_map_.end(); ++pose_inputter_it ) {
			TR << pose_inputter_it->first << ", ";
		}
		TR << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( pose_inputter_type + " is not known to the PoseInputterFactory."
			" Was it registered via a PoseInputterRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return PoseInputterOP();
	}
}

PoseInputSources
PoseInputterFactory::pose_inputs_from_command_line() const
{
	PoseInputSources input_sources;
	for ( PoseInputterMap::const_iterator iter = pose_inputter_creator_map_.begin();
			iter != pose_inputter_creator_map_.end(); ++iter ) {
		// TR << pose_inputter_it->first << ", ";
		PoseInputterOP inputter = iter->second->create_inputter();
		if ( inputter->job_available_on_command_line() ) {
			PoseInputSources iter_sources = inputter->pose_input_sources_from_command_line();
			input_sources.reserve( input_sources.size() + iter_sources.size() );
			for ( core::Size ii = 1; ii <= iter_sources.size(); ++ii ) {
				input_sources.push_back( iter_sources[ ii ] );
			}
		}
	}
	return input_sources;
}

void PoseInputterFactory::set_throw_on_double_registration()
{
	throw_on_double_registration_ = true;
}


void PoseInputterFactory::define_pose_inputter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	utility::tag::define_xml_schema_group(
		pose_inputter_creator_map_,
		pose_inputter_xml_schema_group_name(),
		& PoseInputterFactory::complex_type_name_for_pose_inputter,
		xsd );
}

void PoseInputterFactory::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	for ( CreatorList::const_iterator iter = creator_list_.begin(); iter != creator_list_.end(); ++iter ) {
		(*iter)->list_options_read( read_options );
	}
}


std::string PoseInputterFactory::pose_inputter_xml_schema_group_name()
{
	return "pose_inputter";
}

std::string PoseInputterFactory::complex_type_name_for_pose_inputter( std::string const & inputter_key )
{
	return "pose_inputter_" + inputter_key + "_type";
}

} //namespace pose_inputters
} //namespace jd3
} //namespace protocols

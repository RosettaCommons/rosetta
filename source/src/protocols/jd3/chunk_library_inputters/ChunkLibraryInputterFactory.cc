// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library_inputters/ChunkLibraryInputterFactory.cc
/// @brief  ChunkLibraryInputterFactory class that holds the list of classes able to create Poses from inputs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputterFactory.hh>

#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputter.hh>
#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputterCreator.hh>

// core headers
#include <core/types.hh>

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
namespace chunk_library_inputters {

static basic::Tracer TR( "protocols.jd3.chunk_library_inputters.ChunkLibraryInputterFactory" );

ChunkLibraryInputterFactory::ChunkLibraryInputterFactory() :
	throw_on_double_registration_( false )
{
}

//ChunkLibraryInputterFactory::~ChunkLibraryInputterFactory(){}

/// @brief add a ChunkLibraryInputter prototype, using its default type name as the map key
void
ChunkLibraryInputterFactory::factory_register( ChunkLibraryInputterCreatorOP creator )
{
	runtime_assert( creator != 0 );
	std::string const chunk_library_inputter_type( creator->keyname() );
	if ( chunk_library_inputter_creator_map_.find( chunk_library_inputter_type ) != chunk_library_inputter_creator_map_.end() ) {
		std::string err_msg = "ChunkLibraryInputterFactory::factory_register already has a chunk_library_inputter creator with name \""
			+ chunk_library_inputter_type + "\".  Conflicting ChunkLibraryInputter names";
		if ( throw_on_double_registration_ ) {
			throw CREATE_EXCEPTION( utility::excn::Exception, err_msg );
		} else {
			utility_exit_with_message( err_msg );
		}
	}
	chunk_library_inputter_creator_map_[ chunk_library_inputter_type ] = creator;
	creator_list_.push_back( creator );
}


/// @brief return new ChunkLibraryInputter by key lookup in chunk_library_inputter_prototype_map_ (new ChunkLibraryInputter parses Tag if provided)
ChunkLibraryInputterOP
ChunkLibraryInputterFactory::new_chunk_library_inputter( std::string const & chunk_library_inputter_type ) const
{
	ChunkLibraryInputterMap::const_iterator iter( chunk_library_inputter_creator_map_.find( chunk_library_inputter_type ) );
	if ( iter != chunk_library_inputter_creator_map_.end() ) {
		if ( ! iter->second ) {
			//AMW TODO
			//throw utility::excn::EXCN_RosettaScriptsOption( "Error: ChunkLibraryInputterCreatorOP prototype for " + chunk_library_inputter_type + " is NULL!" );
		}
		return iter->second->create_inputter();
	} else {
		TR<<"Available ChunkLibraryInputters: ";
		for ( auto const & chunk_library_inputter_elem : chunk_library_inputter_creator_map_ ) {
			TR << chunk_library_inputter_elem.first << ", ";
		}
		TR << std::endl;
		// AMW TODO
		// throw utility::excn::EXCN_RosettaScriptsOption( chunk_library_inputter_type + " is not known to the ChunkLibraryInputterFactory."
		// " Was it registered via a ChunkLibraryInputterRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return ChunkLibraryInputterOP();
	}
}

ChunkLibraryInputterFactory::ChunkLibraryInputSourcesAndInputters
ChunkLibraryInputterFactory::chunk_library_inputs_from_command_line() const
{
	ChunkLibraryInputSourcesAndInputters input_sources;
	for ( auto const & elem : chunk_library_inputter_creator_map_ ) {
		ChunkLibraryInputterOP inputter = elem.second->create_inputter();
		if ( inputter->job_available_on_command_line() ) {
			// AMW: there should only be one of these, probably!
			ChunkLibraryInputSources iter_sources = inputter->chunk_library_input_sources_from_command_line();
			runtime_assert( iter_sources.size() == 1 );
			input_sources.reserve( input_sources.size() + iter_sources.size() );
			for ( core::Size ii = 1; ii <= iter_sources.size(); ++ii ) {
				input_sources.push_back( std::make_pair( iter_sources[ ii ], inputter ) );
			}
		}
	}
	return input_sources;
}

void ChunkLibraryInputterFactory::set_throw_on_double_registration()
{
	throw_on_double_registration_ = true;
}


void ChunkLibraryInputterFactory::define_chunk_library_inputter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	utility::tag::define_xml_schema_group(
		chunk_library_inputter_creator_map_,
		chunk_library_inputter_xml_schema_group_name(),
		& ChunkLibraryInputterFactory::complex_type_name_for_chunk_library_inputter,
		xsd );
}

void ChunkLibraryInputterFactory::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	for ( auto const & creator : creator_list_ ) {
		creator->list_options_read( read_options );
	}
}

ChunkLibraryInputterFactory::CreatorList const &
ChunkLibraryInputterFactory::chunk_library_inputter_creators() const
{
	return creator_list_;
}

std::string ChunkLibraryInputterFactory::chunk_library_inputter_xml_schema_group_name()
{
	return "chunk_library_inputter";
}

std::string ChunkLibraryInputterFactory::complex_type_name_for_chunk_library_inputter( std::string const & inputter_key )
{
	return "chunk_library_inputter_" + inputter_key + "_type";
}

} //namespace chunk_library_inputters
} //namespace jd3
} //namespace protocols

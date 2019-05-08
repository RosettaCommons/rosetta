// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/multistage_rosetta_scripts/TagManager.cc
/// @brief
/// @detailed
/// @author Jack Maguire, jackmaguire1444@gmail.com


#include <protocols/multistage_rosetta_scripts/TagManager.hh>

#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/NullMover.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/filters/BasicFilters.hh>

#include <protocols/parser/DataLoaderFactory.hh>
#include <protocols/parser/DataLoader.hh>

#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.multistage_rosetta_scripts.TagManager" );

using namespace utility;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/list.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace multistage_rosetta_scripts {

bool NoFailDataMap::add(
	std::string const & type,
	std::string const & name,
	utility::pointer::ReferenceCountOP const op
) {
	if ( ! has( type, name ) ) { //redundant but avoids TR output that might worry users
		DataMap::add( type, name, op );
	}
	return true;
}

TagManager::TagManager( ) :
	most_recent_request_( 0 )
{}

TagManager::TagManager( core::Size num_input_pose_ids ) :
	most_recent_request_( 0 )
{
	set_num_input_pose_ids( num_input_pose_ids );
}

TagManager::~TagManager(){}


void
TagManager::set_num_input_pose_ids( core::Size num_input_pose_ids ){
	TagListOP default_tag_list( pointer::make_shared< TagList >() );
	tag_list_for_input_pose_id_.assign( num_input_pose_ids + 1, default_tag_list );
}

void
TagManager::register_data_tags_for_input_pose_id(
	core::Size input_pose_id,
	std::list< utility::tag::TagCOP > const & tags_in_order
){
	if ( tags_in_order.size() == 0 ) return;

	TagListOP tag_list( pointer::make_shared< std::list< utility::tag::TagCOP > >( tags_in_order ) );

	auto & common_tags = tag_list_for_input_pose_id_[ 0 ];
	for ( utility::tag::TagCOP tag : * common_tags ) {
		tag_list->push_back( tag );
	}

	tag_list_for_input_pose_id_[ input_pose_id ] = tag_list;
}

void
TagManager::register_data_tags_for_input_pose_id(
	core::Size input_pose_id,
	std::list< utility::tag::TagCOP > && tags_in_order
){
	if ( tags_in_order.size() == 0 ) return;

	TagListOP tag_list( pointer::make_shared< std::list< utility::tag::TagCOP > >( std::move( tags_in_order ) ) );

	auto & common_tags = tag_list_for_input_pose_id_[ 0 ];
	for ( utility::tag::TagCOP tag : * common_tags ) {
		tag_list->push_back( tag );
	}

	tag_list_for_input_pose_id_[ input_pose_id ] = tag_list;
}

ParsedTagCacheOP
TagManager::generate_data_for_input_pose_id(
	core::Size input_pose_id,
	core::pose::Pose const & pose
){
	//We actually want to create a new DataMap every time regardless
	/*if ( most_recent_request_ ) {
	if ( input_pose_id == most_recent_request_->input_pose_id ) {
	return most_recent_request_;
	}
	}*/

	NoFailDataMapOP data_map( pointer::make_shared< NoFailDataMap >() );
	TagMapOP mover_tags = pointer::make_shared< std::map< std::string, utility::tag::TagCOP > >();
	TagMapOP filter_tags = pointer::make_shared< std::map< std::string, utility::tag::TagCOP > >();

	utility::pointer::shared_ptr< moves::Movers_map > movers_map = pointer::make_shared< moves::Movers_map >();
	movers_map->insert( std::pair< std::string, moves::MoverOP >( "null", pointer::make_shared< moves::NullMover >( ) ) );

	utility::pointer::shared_ptr< filters::Filters_map > filters_map = pointer::make_shared< filters::Filters_map >();
	filters::FilterOP true_filter( pointer::make_shared< filters::TrueFilter >() );
	filters::FilterOP false_filter( pointer::make_shared< filters::FalseFilter >() );
	filters_map->insert( std::pair< std::string, filters::FilterOP >( "true_filter", true_filter ) );
	filters_map->insert( std::pair< std::string, filters::FilterOP >( "false_filter", false_filter ) );

	moves::MoverFactory const * mover_factory = moves::MoverFactory::get_instance();
	filters::FilterFactory const * filter_factory = filters::FilterFactory::get_instance();

	for ( utility::tag::TagCOP tag : * tag_list_for_input_pose_id_[ input_pose_id ] ) {
		if ( tag->getName() == "MOVERS" ) {
			for ( auto & mover_tag : tag->getTags() ) {
				std::string const mover_name = mover_tag->getOption< std::string >( "name" );
				( * mover_tags )[ mover_name ] = mover_tag;

				movers_map->insert(
					std::pair< std::string, moves::MoverOP >( mover_name,
					mover_factory->newMover( mover_tag, * data_map, * filters_map, * movers_map, pose )
					)
				);
			}
		} else if ( tag->getName() == "FILTERS" ) {
			for ( auto & filter_tag : tag->getTags() ) {
				std::string const filter_name = filter_tag->getOption< std::string >( "name" );
				( * filter_tags )[ filter_name ] = filter_tag;

				filters_map->insert(
					std::pair< std::string, filters::FilterOP >( filter_name,
					filter_factory->newFilter( filter_tag, * data_map, * filters_map, * movers_map, pose )
					)
				);
			}

		} else {
			parser::DataLoaderOP loader = parser::DataLoaderFactory::get_instance()->newDataLoader( tag->getName() );
			loader->load_data( pose, tag, * data_map );
		}
	}

	most_recent_request_ = pointer::make_shared< ParsedTagCache >( input_pose_id, data_map, mover_tags, filter_tags, movers_map, filters_map );
	return most_recent_request_;
}


} //multistage_rosetta_scripts
} //protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::TagManager::save( Archive & arc ) const {
	arc( CEREAL_NVP( tag_list_for_input_pose_id_ ) ); // std::vector<TagListOP>
	//arc( CEREAL_NVP( most_recent_request_ ) ); // ParsedTagCacheOP
	// EXEMPT most_recent_request_
	arc( CEREAL_NVP( dummy_pose_ ) ); // core::pose::Pose
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::TagManager::load( Archive & arc ) {

	std::vector< std::shared_ptr< std::list< utility::tag::TagOP > > > local_tag_list_for_input_pose_id;
	arc( local_tag_list_for_input_pose_id ); // std::vector<TagListOP>
	std::vector< std::shared_ptr< std::list< utility::tag::TagCOP > > > tmp;
	tmp.reserve(local_tag_list_for_input_pose_id.size());
	for ( auto taglist: local_tag_list_for_input_pose_id ) {
		auto tags = std::make_shared< std::list< utility::tag::TagCOP > >();
		for ( auto tag: *taglist ) {
			tags->push_back(tag);
		}
		tmp.push_back(tags);
	}
	tag_list_for_input_pose_id_ = tmp; // copy the non-const pointer(s) into the const pointer(s)

	//arc( most_recent_request_ ); // ParsedTagCacheOP
	// EXEMPT most_recent_request_
	arc( dummy_pose_ ); // core::pose::Pose
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::multistage_rosetta_scripts::TagManager );
CEREAL_REGISTER_TYPE( protocols::multistage_rosetta_scripts::TagManager )


/////// @brief Automatically generated serialization method
////template< class Archive >
////void
////protocols::multistage_rosetta_scripts::NoFailDataMap::save( Archive & arc ) const {
//// arc( cereal::base_class< basic::datacache::DataMap >( this ) );
////}
////
/////// @brief Automatically generated deserialization method
////template< class Archive >
////void
////protocols::multistage_rosetta_scripts::NoFailDataMap::load( Archive & arc ) {
//// arc( cereal::base_class< basic::datacache::DataMap >( this ) );
////}
////
////SAVE_AND_LOAD_SERIALIZABLE( protocols::multistage_rosetta_scripts::NoFailDataMap );
////CEREAL_REGISTER_TYPE( protocols::multistage_rosetta_scripts::NoFailDataMap )


/// @brief Automatically generated serialization method
//// template< class Archive >
//// void
//// protocols::multistage_rosetta_scripts::ParsedTagCache::save( Archive & arc ) const {
////  arc( CEREAL_NVP( input_pose_id ) ); // core::Size
////  arc( CEREAL_NVP( data_map ) ); // NoFailDataMapOP
////  arc( CEREAL_NVP( mover_tags ) ); // TagMapOP
////  arc( CEREAL_NVP( filter_tags ) ); // TagMapOP
////  arc( CEREAL_NVP( movers_map ) ); // utility::pointer::shared_ptr<moves::Movers_map>
////  arc( CEREAL_NVP( filters_map ) ); // utility::pointer::shared_ptr<filters::Filters_map>
//// }
////
//// /// @brief Automatically generated deserialization method
//// template< class Archive >
//// void
//// protocols::multistage_rosetta_scripts::ParsedTagCache::load( Archive & arc ) {
////  arc( input_pose_id ); // core::Size
////  arc( data_map ); // NoFailDataMapOP
////  std::shared_ptr< std::map< std::basic_string< char >, std::shared_ptr< utility::tag::Tag >, struct std::less< std::basic_string< char > > > > local_mover_tags;
////  arc( local_mover_tags ); // TagMapOP
////  mover_tags = local_mover_tags; // copy the non-const pointer(s) into the const pointer(s)
////  std::shared_ptr< std::map< std::basic_string< char >, std::shared_ptr< utility::tag::Tag >, struct std::less< std::basic_string< char > > > > local_filter_tags;
////  arc( local_filter_tags ); // TagMapOP
////  filter_tags = local_filter_tags; // copy the non-const pointer(s) into the const pointer(s)
////  std::shared_ptr< std::map< std::basic_string< char >, std::shared_ptr< protocols::moves::Mover >, struct std::less< std::basic_string< char > > > > local_movers_map;
////  arc( local_movers_map ); // utility::pointer::shared_ptr<moves::Movers_map>
////  movers_map = local_movers_map; // copy the non-const pointer(s) into the const pointer(s)
////  std::shared_ptr< std::map< std::basic_string< char >, std::shared_ptr< protocols::filters::Filter >, struct std::less< std::basic_string< char > > > > local_filters_map;
////  arc( local_filters_map ); // utility::pointer::shared_ptr<filters::Filters_map>
////  filters_map = local_filters_map; // copy the non-const pointer(s) into the const pointer(s)
//// }
////
//// SAVE_AND_LOAD_SERIALIZABLE( protocols::multistage_rosetta_scripts::ParsedTagCache );
CEREAL_REGISTER_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_TagManager )
#endif // SERIALIZATION

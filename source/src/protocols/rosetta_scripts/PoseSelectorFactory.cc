// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Singleton instance and mutex static data members
namespace utility {

using protocols::rosetta_scripts::PoseSelectorFactory;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< PoseSelectorFactory >::singleton_mutex_{};
template <> std::atomic< PoseSelectorFactory * > utility::SingletonBase< PoseSelectorFactory >::instance_( 0 );
#else
template <> PoseSelectorFactory * utility::SingletonBase< PoseSelectorFactory >::instance_( 0 );
#endif

}

namespace protocols {
namespace rosetta_scripts {

static thread_local basic::Tracer TR( "protocols.rosetta_scripts.PoseSelectorFactory" );

PoseSelectorFactory *
PoseSelectorFactory::create_singleton_instance()
{
	return new PoseSelectorFactory;
}

PoseSelectorFactory::PoseSelectorFactory(){}

PoseSelectorFactory::~PoseSelectorFactory(){}

/// @brief add a PoseSelector prototype, using its default type name as the map key
void
PoseSelectorFactory::factory_register( PoseSelectorCreatorOP creator )
{
	runtime_assert( creator != 0 );
	std::string const pose_selector_type( creator->keyname() );
	if ( pose_selector_type == "UNDEFINED NAME" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Can't map derived PoseSelector with undefined type name.");
	}
	if ( poseselector_creator_map_.find( pose_selector_type ) != poseselector_creator_map_.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption("PoseSelectorFactory::factory_register already has a pose selector creator with name \"" + pose_selector_type + "\".  Conflicting pose selector names" );
	}
	poseselector_creator_map_[ pose_selector_type ] = creator;
}


/// @brief return new PoseSelector by key lookup in poseselector_creator_map_ (new PoseSelector parses Tag if provided)
PoseSelectorOP
PoseSelectorFactory::newPoseSelector(	std::string const & pose_selector_type )
{
	PoseSelectorMap::const_iterator iter( poseselector_creator_map_.find( pose_selector_type ) );
	if ( iter != poseselector_creator_map_.end() ) {
		if ( ! iter->second ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: PoseSelectorCreatorOP prototype for " + pose_selector_type + " is NULL!" );
		}
		return iter->second->create_selector();
	}
	else {
		TR<<"Available pose selectors: ";
		for( PoseSelectorMap::const_iterator it = poseselector_creator_map_.begin(); it != poseselector_creator_map_.end(); ++it )
			TR<<it->first<<", ";
		TR<<std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( pose_selector_type + " is not known to the PoseSelectorFactory. Was it registered via a PoseSelectorRegistrator in one of the init.cc files?" );
		return NULL;
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
	runtime_assert( selector != 0 );
	selector->parse_my_tag( tag, data, filters, movers, pose );
	return selector;
}

} //namespace rosetta_scripts
} //namespace protocols

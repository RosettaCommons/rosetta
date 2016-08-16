// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/GridFactory.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>
//#include <protocols/qsar/qsarTypeManager.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Singleton instance and mutex static data members
namespace utility {

using protocols::qsar::scoring_grid::GridFactory;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< GridFactory >::singleton_mutex_{};
template <> std::atomic< GridFactory * > utility::SingletonBase< GridFactory >::instance_( 0 );
#else
template <> GridFactory * utility::SingletonBase< GridFactory >::instance_( 0 );
#endif

}

namespace protocols {
namespace qsar {
namespace scoring_grid {

GridFactory *
GridFactory::create_singleton_instance()
{
	return new GridFactory;
}

GridFactory::GridFactory()
{
	//This isn't strictly necessary but I had an absolutely amazing bug involving
	//The gridfactory being prematurely destroyed and then immediately recreated
	//This was resulting in data corruption and it would have been 1000x easier to debug
	//if I had this clear statement here, so I'm leaving it in for the next time.
	grid_creator_map_.clear();
}

GridFactory::~GridFactory()
{ }


/// @brief add a Grid prototype, using it's default type name as the map key
void
GridFactory::factory_register(GridCreatorOP creator)
{
	runtime_assert(creator != 0 );
	std::string grid_type( creator->keyname());
	//qsarType grid_enum = qsarTypeManager::qsar_type_from_name(grid_type);
	if ( grid_type == "UNDEFINED NAME" ) {
		utility_exit_with_message("can't map derived grid with undefined type name.");
	}
	if ( grid_creator_map_.find(grid_type) != grid_creator_map_.end() ) {
		utility_exit_with_message("GridFactory::factory_register already has a move creator with name \"" + grid_type + "\". Conflicting grid names");
	}
	grid_creator_map_.insert(std::pair<std::string,GridCreatorOP>(grid_type,creator));
	grid_creator_map_[grid_type]  = creator;
}

GridBaseOP GridFactory::new_grid(utility::tag::TagCOP tag) const
{
	//string const name = tag->getName();
	std::string const type = tag->getOption<std::string>("grid_type");

	GridMap::const_iterator iter(grid_creator_map_.find(type));
	if ( iter != grid_creator_map_.end() ) {
		if ( !iter->second ) {
			utility_exit_with_message("Error: GridCreatorOP prototype for "+type+ " is NULL!");
			return NULL;
		}
		return iter->second->create_grid(tag);
	} else {
		utility_exit_with_message(type + " is not known to the GridFactory.  Was it registered via a GridRegistrator in one of the init.cc files");
		return NULL;
	}

}

GridBaseOP GridFactory::new_grid(utility::json_spirit::mObject data ) const
{
	//figure out what kind of grid to make.  There will be a "type" tag either in the top level of the
	//heirarchy or in the base_data level.

	std::string type;

	utility::json_spirit::mObject::iterator type_it(data.find("type"));
	if ( type_it != data.end() ) {
		//If this is a metagrid then we can find a type tag in the top level
		type = type_it->second.get_str();

	} else {
		//OK, it's not a metagrid, Everything that inherits from SingleGrid has a "base_data" tag, and "type" is under that.
		utility::json_spirit::mObject base_data(data["base_data"].get_obj());
		type = base_data["type"].get_str();
	}

	//make a new grid
	GridBaseOP new_grid;
	GridMap::const_iterator iter(grid_creator_map_.find(type));
	if ( iter != grid_creator_map_.end() ) {
		if ( !iter->second ) {
			utility_exit_with_message("Error: GridCreatorOP prototype for "+type+ " is NULL!");
		}
		new_grid = iter->second->create_grid();
	} else {
		utility_exit_with_message(type + " is not known to the GridFactory.  Was it registered via a GridRegistrator in one of the init.cc files");
	}

	//deserialize object into new grid and return
	assert(new_grid);
	new_grid->deserialize(data);
	return new_grid;
}

}
}
}

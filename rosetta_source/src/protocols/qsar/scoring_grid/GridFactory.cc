// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /git/src/protocols/qsar/scoring_grid/GridFactory.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
//#include <protocols/qsar/qsarTypeManager.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {

GridFactory * GridFactory::instance_(0);

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

GridFactory *
GridFactory::get_instance()
{
	if(!instance_)
	{
		GridFactory * instance_local = new GridFactory;
		instance_ = instance_local;
	}
	return instance_;
}


///@brief add a Grid prototype, using it's default type name as the map key
void
GridFactory::factory_register(GridCreatorOP creator)
{
	runtime_assert(creator);
	std::string grid_type( creator->keyname());
	//qsarType grid_enum = qsarTypeManager::qsar_type_from_name(grid_type);
	if(grid_type == "UNDEFINED NAME")
	{
		utility_exit_with_message("can't map derived grid with undefined type name.");
	}
	if (grid_creator_map_.find(grid_type) != grid_creator_map_.end())
	{
		utility_exit_with_message("GridFactory::factory_register already has a move creator with name \"" + grid_type + "\". Conflicting grid names");
	}
	grid_creator_map_.insert(std::pair<std::string,GridCreatorOP>(grid_type,creator));
	grid_creator_map_[grid_type]  = creator;
}

GridBaseOP GridFactory::new_grid(utility::tag::TagPtr const tag) const
{
	//string const name = tag->getName();
	std::string const type = tag->getOption<std::string>("grid_type");

	GridMap::const_iterator iter(grid_creator_map_.find(type));
	if( iter != grid_creator_map_.end())
	{
		if(!iter->second)
		{
			utility_exit_with_message("Error: GridCreatorOP prototype for "+type+ " is NULL!");
			return NULL;
		}
		return iter->second->create_grid(tag);
	}
	else
	{
		utility_exit_with_message(type + " is not known to the GridFactory.  Was it registered via a GridRegistrator in one of the init.cc files");
		return NULL;
	}

}

}
}
}

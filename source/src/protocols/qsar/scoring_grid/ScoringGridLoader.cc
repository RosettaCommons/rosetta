// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/qsar/scoring_grid/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/qsar/scoring_grid/ScoringGridLoader.hh>
#include <protocols/qsar/scoring_grid/ScoringGridLoaderCreator.hh>

// Project Headers
#include <protocols/qsar/scoring_grid/GridManager.hh>
//#include <protocols/qsar/qsarTypeManager.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#define foreach BOOST_FOREACH

namespace protocols {
namespace qsar {
namespace scoring_grid {

static basic::Tracer TR( "protocols.qsar.scoring_grid.ScoringGridLoader" );

ScoringGridLoader::ScoringGridLoader() {}
ScoringGridLoader::~ScoringGridLoader() {}

void ScoringGridLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap &
) const
{
	using namespace utility::tag;
	typedef utility::vector0< TagCOP > TagCOPs;

	/// Setup the scoring grid_manager

	//core::Real width = 40.0;
	//core::Real resolution = 0.25;

	qsar::scoring_grid::GridManager* grid_manager(qsar::scoring_grid::GridManager::get_instance());

	if( tag->hasOption("width") ) {
		grid_manager->set_width(tag->getOption<core::Real>("width"));

	}
	if( tag->hasOption("resolution") ){
		grid_manager->set_resolution(tag->getOption<core::Real>("resolution"));
	}

	if(tag->hasOption("ligand_chain") ) {
		grid_manager->set_chain(tag->getOption<char>("ligand_chain"));
	}

	if(tag->hasOption("normalize_mode"))
	{
		std::string normalize_mode = tag->getOption<std::string>("normalize_mode");
		grid_manager->set_normalization_function(normalize_mode);
	}

	/// Add grids to the scoring grid manager

	TagCOPs const grid_tags( tag->getTags() );
	if (grid_tags.size()==0){
		TR <<"WARNING WARNING grid manager will be empty" <<std::endl;
	}

	foreach(TagCOP tag, grid_tags){
		grid_manager->make_new_grid(tag);
	}

	TR.flush();
}

protocols::jd2::parser::DataLoaderOP
ScoringGridLoaderCreator::create_loader() const { return new ScoringGridLoader; }

std::string
ScoringGridLoaderCreator::keyname() const { return "SCORINGGRIDS"; }


} //namespace scoring_grid
} //namespace qsar
} //namespace protocols

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/jd2/parser/ScoringGridLoader.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>

// Project Headers
#include <protocols/qsar/scoring_grid/GridManager.hh>
//#include <protocols/qsar/qsarTypeManager.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>


namespace protocols {
namespace jd2 {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.ScoringGridLoader" );

ScoringGridLoader::ScoringGridLoader() {}
ScoringGridLoader::~ScoringGridLoader() {}

void ScoringGridLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagPtr const tag,
	moves::DataMap & data
) const
{
	using namespace utility::tag;
	typedef utility::vector0< TagPtr > TagPtrs;

	TagPtrs const TO_tags( tag->getTags() );

	core::Real width = 40.0;
	core::Real resolution = 0.25;
	utility::vector1<std::string> grids_to_add;

	for( TagPtrs::const_iterator tp( TO_tags.begin() ), tp_e( TO_tags.end() ); tp != tp_e; ++tp ) {
		std::string const type( (*tp)->getName());
		if(type == "dimensions") {
			if((*tp)->hasOption("width") && (*tp)->hasOption("resolution")) {
				width = (*tp)->getOption<core::Real>("width");
				resolution = (*tp)->getOption<core::Real>("resolution");
			} else {
				utility_exit_with_message("you must specify both width and resolution");
			}
		} else if(type=="grid") {
			if ( ! (*tp)->hasOption("grid_type") ) {
				utility_exit_with_message("you need to specify a grid type");
			}
			std::string grid_name((*tp)->getOption<std::string>("grid_type"));
			grids_to_add.push_back(grid_name);
		}
	}

	if(grids_to_add.size() == 0) {
		TR <<"WARNING WARNING grid manager will be empty" <<std::endl;
	}

	qsar::scoring_grid::GridManagerOP grid_manager(qsar::scoring_grid::GridManager::get_instance());
	grid_manager->set_dimensions(width, resolution);
	//qsar::scoring_grid::GridManagerOP grid_manager(new qsar::scoring_grid::GridManager(width,resolution));
	utility::vector1<std::string>::iterator current_grid;
	for(current_grid = grids_to_add.begin();current_grid != grids_to_add.end(); ++current_grid) {
		grid_manager->make_new_grid(*current_grid);
		TR.Debug <<"adding grid: " <<*current_grid <<std::endl;
	}

	data.add( "scoringgrid", "default", grid_manager );
	TR.flush();
}

DataLoaderOP
ScoringGridLoaderCreator::create_loader() const { return new ScoringGridLoader; }

std::string
ScoringGridLoaderCreator::keyname() const { return "SCORINGGRIDS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols

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

	/// Add grids to the scoring grid manager

	TagPtrs const grid_tags( tag->getTags() );
	if (grid_tags.size()==0){
		TR <<"WARNING WARNING grid manager will be empty" <<std::endl;
	}

	for( TagPtrs::const_iterator tp( grid_tags.begin() ), tp_e( grid_tags.end() ); tp != tp_e; ++tp ) {
		grid_manager->make_new_grid(*tp);
	}

	TR.flush();
}

DataLoaderOP
ScoringGridLoaderCreator::create_loader() const { return new ScoringGridLoader; }

std::string
ScoringGridLoaderCreator::keyname() const { return "SCORINGGRIDS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ScoreMap.cc
/// @brief  A place to put some common functions for scoremap output
/// @author Monica Berrondo

// Unit headers
#include <protocols/jd2/ScoreMap.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>


// ObjexxFCL headers

#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


using namespace core;

namespace protocols {
namespace jd2 {

/// @details Auto-generated virtual destructor
ScoreMap::~ScoreMap() {}

/// Score output helper functions
/// @details
///  creates a list of non zero weighted energies and adds them to the
///  score_file information that is to be written out at the end of the
///  protocol.
void
ScoreMap::nonzero_energies(
	std::map < std::string, core::Real > & score_map,
	scoring::ScoreFunctionOP score_fxn,
	pose::Pose & pose
)
{
	using namespace core::scoring;

	(*score_fxn)(pose);

	score_map_from_scored_pose(score_map, pose);
}

/// @details creates score map from scored pdb; const so it can be used in job distributor
void
ScoreMap::score_map_from_scored_pose(
	std::map < std::string, core::Real > & score_map,
	pose::Pose const & pose
) {
	using namespace core::scoring;

	// Which score terms to use
	core::scoring::EnergyMap weights = pose.energies().weights();
	typedef utility::vector1<core::scoring::ScoreType> ScoreTypeVec;
	ScoreTypeVec score_types;
	for ( int i = 1; i <= core::scoring::n_score_types; ++i ) {
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);
		if ( weights[ii] != 0 ) score_types.push_back(ii);
	}

	core::Real total(0);

	for ( ScoreTypeVec::iterator ii = score_types.begin(), end_ii = score_types.end(); ii != end_ii; ++ii ) {
		core::Real const some_score( weights[(*ii)] * pose.energies().total_energies()[ *ii ] );
		score_map[ name_from_score_type(*ii) ] = some_score;
		total += some_score;
	}
	score_map[ name_from_score_type(core::scoring::total_score) ] = total;
}

std::map< std::string, core::Real > ScoreMap::score_map_from_scored_pose( core::pose::Pose const & pose ){
	std::map< std::string, core::Real > score_map;
	score_map_from_scored_pose(score_map, pose);
	return score_map;
}


/// @brief print out the contents of the ScoreMap
void
ScoreMap::print(
	std::map < std::string, core::Real > & score_map,
	std::ostream & out
)
{
	using namespace ObjexxFCL::format;

	std::map< std::string, Real >::const_iterator pair;
	Size width (8), precision (3);

	//print the header
	for ( pair=score_map.begin(); pair!=score_map.end(); ++pair ) {
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << A( width, pair->first );
	}
	out << std::endl;

	//print the information
	for ( pair=score_map.begin(); pair!=score_map.end(); ++pair ) {
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << F( width, precision, pair->second );
	}
	out << std::endl;

}


} //jd2
} //protocols

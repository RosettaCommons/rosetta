// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ScoreMap.hh
/// @brief  A place to put some common functions for scoremap output
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_jd2_ScoreMap_hh
#define INCLUDED_protocols_jd2_ScoreMap_hh

#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <map>


namespace protocols {
namespace jd2 {

class ScoreMap : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ScoreMap() override;
	/// @brief full atom energies for output
	static void nonzero_energies(
		std::map< std::string, core::Real> & score_map,
		core::scoring::ScoreFunctionOP score_fxn,
		core::pose::Pose & pose
	);

	/// @brief generates a scoremap assuming the pose is already scored (note const w.r.t. pose)
	static void score_map_from_scored_pose(
		std::map< std::string, core::Real> & score_map,
		core::pose::Pose const & pose
	);

	/// @brief return-by-value version of score_map_from_scored_pose
	static std::map< std::string, core::Real> score_map_from_scored_pose( core::pose::Pose const & pose );

	/// @brief print out the values in the scoremap
	static void print(
		std::map < std::string, core::Real > & score_map,
		std::ostream & out
	);
};
}// jd2
}// protocols

#endif

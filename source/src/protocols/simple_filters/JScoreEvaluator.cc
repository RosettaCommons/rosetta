// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/JScoreEvaluator.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/util/SwitchResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_filters/JScoreEvaluator.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols  {
namespace simple_filters {

JScoreEvaluator::JScoreEvaluator(
	std::string const & weights,
	std::string const & type_set_name
) :
	evaluation::SingleValuePoseEvaluator< core::Real >( weights + "_" + type_set_name ),
	scorefxn_(core::scoring::ScoreFunctionFactory::create_score_function(weights)),
	type_set_name_(type_set_name),
	col_name_(weights)
{}

JScoreEvaluator::~JScoreEvaluator() {}

void JScoreEvaluator::apply(
	core::pose::Pose & pose_in,
	std::string /*tag*/,
	core::io::silent::SilentStruct & ss
) const {
	// make a copy of the Pose to prevent nuking side-chains when switching
	// ResidueTypes
	core::pose::Pose pose(pose_in);

	// switch residue types
	core::util::switch_to_residue_type_set( pose, type_set_name_);

	// score with ScoreFunction
	core::Real const score( (*scorefxn_)(pose) );

	// put score into SilentStruct
	ss.add_energy( col_name_, score );
}

} // simple_filter
} // protocols

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/evaluation/JScoreEvaluator.hh
/// @brief
/// @author James Thompson
///

#ifndef _INCLUDED_protocols_evaluation_JScoreEvaluator_hh_
#define _INCLUDED_protocols_evaluation_JScoreEvaluator_hh_

#include <string>
#include <utility/vector1.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

namespace protocols {
namespace evaluation {

class JScoreEvaluator : public SingleValuePoseEvaluator< core::Real > {
public:
	JScoreEvaluator(
		std::string const & weights       = "score12_full",
		std::string const & type_set_name = "fa_standard"
	);

	~JScoreEvaluator();

	virtual void apply(
		core::pose::Pose & pose,
		std::string tag,
		core::io::silent::SilentStruct & ss
	) const;

	virtual core::Real apply(
		core::pose::Pose & /*pose*/
	) const {
		return 0;
	}

private:
	 core::scoring::ScoreFunctionOP scorefxn_;
	 std::string const type_set_name_;
	 std::string const col_name_;
};

} // evaluation
} // protocols

#endif  // _INCLUDED_protocols_evaluation_JScoreEvaluator_hh_

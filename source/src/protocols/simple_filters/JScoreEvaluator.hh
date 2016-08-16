// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/JScoreEvaluator.hh
/// @brief
/// @author James Thompson


#ifndef INCLUDED_protocols_simple_filters_JScoreEvaluator_hh
#define INCLUDED_protocols_simple_filters_JScoreEvaluator_hh

#include <string>
#include <utility/vector1.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

namespace protocols {
namespace simple_filters {

class JScoreEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
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

} // simple_filter
} // protocols

#endif  // _INCLUDED_protocols_simple_filters_JScoreEvaluator_hh_

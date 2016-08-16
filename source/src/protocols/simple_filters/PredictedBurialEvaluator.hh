// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/BurialEvaluation.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_simple_filters_PredictedBurialEvaluator_hh
#define INCLUDED_protocols_simple_filters_PredictedBurialEvaluator_hh


#include <string>
#include <utility/vector1.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

namespace protocols {
namespace simple_filters {

class PredictedBurialEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	PredictedBurialEvaluator(
		std::string const & fn
	);

	~PredictedBurialEvaluator();

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
	void init_from_file( std::string const & fn );
	utility::vector1< core::Real > predicted_burial_;
};

} // simple_filter
} // protocols

#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/StructuralSimilarityEvaluation.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_simple_filters_StructuralSimilarityEvaluator_hh
#define INCLUDED_protocols_simple_filters_StructuralSimilarityEvaluator_hh


#include <string>
#include <utility/vector1.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

namespace protocols {
namespace simple_filters {

class StructuralSimilarityEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	StructuralSimilarityEvaluator(
		utility::vector1< core::pose::Pose > const & poses,
		std::string const & atom_name = "CA",
		std::string const & tag = "sim"
	);

	~StructuralSimilarityEvaluator();

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
	 // KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	 // std::string const & atom_name_;
	 utility::vector1< core::pose::Pose > const poses_;
};

} // simple_filter
} // protocols

#endif

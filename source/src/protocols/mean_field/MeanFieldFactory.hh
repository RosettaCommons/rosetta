// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/mean_field/MeanFieldFactory.hh
/// @brief  Mean field calculator class declaration
/// @author Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_MeanFieldFactory_hh
#define INCLUDED_protocols_mean_field_MeanFieldFactory_hh

// Unit headers
#include <protocols/mean_field/MeanFieldFactory.fwd.hh>

// Package headers
#include <protocols/mean_field/MeanField.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <utility/vector1.hh>
#include <core/types.hh>

namespace protocols {
namespace mean_field {

class MeanFieldFactory {
public:

	/// @brief initializes MeanField class based on options
	static
	MeanFieldOP
	create_mean_field(
		core::Size option,
		core::pose::PoseOPs & poses,
		utility::vector1 < core::pack::task::PackerTaskOP > tasks,
		core::scoring::ScoreFunctionOP scfxn,
		core::Real lambda_mem,
		core::Real conv_diff,
		core::Real temp,
		core::Real threshold
	);
};

}
}

#endif

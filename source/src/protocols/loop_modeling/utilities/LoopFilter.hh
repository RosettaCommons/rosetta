// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_LoopFilter_HH
#define INCLUDED_protocols_loop_modeling_utilities_LoopFilter_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/LoopFilter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Apply any filter in the context of a loop modeling simulation.
class LoopFilter : public LoopMover {

public:
	/// @brief Constructor with filter argument.
	LoopFilter(protocols::filters::FilterOP filter);

	/// @copydoc LoopMover::get_name
	string get_name() const { return "LoopFilter"; }

public:
	/// @brief Apply the given filter and return the result.
	bool do_apply(Pose & pose);

private:
	protocols::filters::FilterOP filter_;

};

}
}
}


#endif


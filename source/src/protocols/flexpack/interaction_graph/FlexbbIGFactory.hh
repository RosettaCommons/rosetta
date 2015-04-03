// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/flexpack/interaction_graph/FlexbbIGFactory.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Florian Richter (floric@u.washington.edu)

#ifndef INCLUDED_protocols_flexpack_interaction_graph_FlexbbIGFactory_hh
#define INCLUDED_protocols_flexpack_interaction_graph_FlexbbIGFactory_hh

#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.fwd.hh>

// Package headers
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.fwd.hh>

// Project headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {
namespace interaction_graph {

class FlexbbIGFactory{

public:

	static
	FlexbbInteractionGraphOP
	create_flexbb_interaction_graph(
		core::pack::task::PackerTask const & task,
		rotamer_set::FlexbbRotamerSets const & flexrots,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn
	);

};

}
}
}

#endif

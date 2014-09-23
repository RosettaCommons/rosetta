// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/flexpack/
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <protocols/flexpack/interaction_graph/FlexbbIGFactory.hh>

/// Package headers
#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.hh>
#include <protocols/flexpack/interaction_graph/MinimalistFlexbbInteractionGraph.hh>
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {
namespace interaction_graph{


FlexbbInteractionGraphOP
FlexbbIGFactory::create_flexbb_interaction_graph(
	core::pack::task::PackerTask const & /*task*/,
	rotamer_set::FlexbbRotamerSets const & flexrots,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & scorefxn
)
{
	MinimalistFlexbbInteractionGraphOP ig( new MinimalistFlexbbInteractionGraph( flexrots.nmoltenres() ) );
	ig->set_pose( pose );
	ig->set_scorefxn( scorefxn );
	return ig;

}


}
}
}



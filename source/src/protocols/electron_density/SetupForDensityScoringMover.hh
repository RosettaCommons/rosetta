// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @detailed
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_electron_density_SetupForDensityScoringMover_hh
#define INCLUDED_protocols_electron_density_SetupForDensityScoringMover_hh

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/electron_density/util.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/loops/Loops.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace electron_density {

// set a pose for density scoring
// wraps 'addVirtualResAsRoot' + 'dockPoseIntoMap'
class SetupForDensityScoringMover : public moves::Mover {
public:
	SetupForDensityScoringMover();

	virtual void apply( core::pose::Pose & pose );

	moves::MoverOP clone() const;

	virtual std::string get_name() const;

	virtual void mask( protocols::loops::Loops const & loops );

	core::Real getScore() { return last_score; }

	virtual void parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap &data,
			filters::Filters_map const &filters,
			moves::Movers_map const &movers,
			core::pose::Pose const & pose );

private:
	std::string dock_into_dens_strategy_;
	utility::vector1< core::Size > mask_reses_;
	core::Real last_score;
};

// set up fold tree _and_ scoring function for density scoring
void set_pose_and_scorefxn_for_edens_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction &scorefxn );

typedef utility::pointer::owning_ptr< SetupForDensityScoringMover > SetupForDensityScoringMoverOP;

}
}

#endif

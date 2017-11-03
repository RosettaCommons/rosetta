// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_electron_density_SetupForDensityScoringMover_hh
#define INCLUDED_protocols_electron_density_SetupForDensityScoringMover_hh

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

//// C++ headers
#include <string>

#include <protocols/loops/Loops.fwd.hh>

#include <utility/vector1.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

namespace protocols {
namespace electron_density {

// set a pose for density scoring
// wraps 'addVirtualResAsRoot' + 'dockPoseIntoMap'
class SetupForDensityScoringMover : public moves::Mover {
public:
	SetupForDensityScoringMover();
	SetupForDensityScoringMover( utility::options::OptionCollection const & options );

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;


	virtual void mask( protocols::loops::Loops const & loops );

	core::Real getScore() { return last_score; }

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static
	void
	options_read_in_ctor( utility::options::OptionKeyList & opts );

private:
	std::string dock_into_dens_strategy_;
	utility::vector1< core::Size > mask_reses_;
	core::Real last_score;
};

// set up fold tree _and_ scoring function for density scoring
void set_pose_and_scorefxn_for_edens_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction &scorefxn );

typedef utility::pointer::shared_ptr< SetupForDensityScoringMover > SetupForDensityScoringMoverOP;

}
}

#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CartminLoopCloser.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_devel_loop_creation_CartminLoopCloser_HH
#define INCLUDED_devel_loop_creation_CartminLoopCloser_HH

//Unit
#include <devel/loop_creation/LoopCloser.hh>

//core
#include <core/scoring/ScoreFunction.fwd.hh>

//protocols
#include <protocols/loops/Loop.hh>

namespace devel {
namespace loop_creation {

class CartminLoopCloser : public LoopCloser
{
public:

	/// @brief default constructor
	CartminLoopCloser();

	/// @brief explicit constructor
	CartminLoopCloser(
		core::scoring::ScoreFunctionOP scorefxn,
		core::Real minimization_tolerance,
		core::Real max_chainbreak
	);

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;

	std::string
	get_name() const;

	void
	init();

	void
	apply ( core::pose::Pose & pose );

	bool
	check_closure ( core::pose::Pose & pose );

	/// @brief parse tag for use in RosettaScripts
	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
	);

private:

	//Should we set up a fold tree that prevents modification outside the loop region?
	bool prevent_nonloop_modifications_;

	core::scoring::ScoreFunctionOP scorefxn_;

	//The threshold for whether or not we consider this loop closed
	core::Real minimization_tolerance_;
	core::Real max_chainbreak_;

};

} //loop creation
} //devel

#endif

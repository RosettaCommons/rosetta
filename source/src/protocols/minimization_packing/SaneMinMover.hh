// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/minimization_packing/SaneMinMover.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_minimization_packing_SaneMinMover_hh
#define INCLUDED_protocols_minimization_packing_SaneMinMover_hh

#include <protocols/minimization_packing/SaneMinMover.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace minimization_packing {

///////////////////////////////////////////////////////////////////////////
// @brief A protocols::moves::Mover that minimizes a Pose to a local energy minimum by
/// performing energy minimization of a ScoreFunction over the allowable degrees
/// of freedom, defined by a MoveMap. Unlike the classic MinMover, the only
/// method for setting Minimization options is via the MinimizerOptions class.
class SaneMinMover : public protocols::moves::Mover {

public:
	SaneMinMover();
	SaneMinMover(std::string const & name);
	SaneMinMover(
		core::kinematics::MoveMapOP movemap_in,
		core::scoring::ScoreFunctionOP scorefxn_in,
		core::optimization::MinimizerOptionsOP min_options_in,
		bool cartesian = false
	);

	~SaneMinMover() override;
	protocols::moves::MoverOP clone() const override;

	// getters
	bool cartesian() const;
	core::kinematics::MoveMapOP move_map() const;
	core::scoring::ScoreFunctionOP score_function() const;
	core::optimization::MinimizerOptionsOP min_options() const;

	/// @brief Minimizes the DOFs of pose specified in the MoveMap
	void apply( core::pose::Pose & pose ) override;


	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	// set reasonable defaults for scorefxn_, movemap_ and min_options_
	void set_defaults_();

	bool cartesian_;
	core::kinematics::MoveMapOP movemap_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::optimization::MinimizerOptionsOP min_options_;
};

} // moves
} // protocols

#endif

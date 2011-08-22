// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/SaneMinMover.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_moves_SaneMinMover_hh
#define INCLUDED_protocols_moves_SaneMinMover_hh

#include <protocols/moves/SaneMinMover.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace moves {

///////////////////////////////////////////////////////////////////////////
// @brief A Mover that minimizes a Pose to a local energy minimum by
/// performing energy minimization of a ScoreFunction over the allowable degrees
/// of freedom, defined by a MoveMap. Unlike the classic MinMover, the only
/// method for setting Minimization options is via the MinimizerOptions class.
class SaneMinMover : public Mover {

	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::optimization::MinimizerOptionsOP MinimizerOptionsOP;
	typedef core::Real Real;

public:
	SaneMinMover();
	SaneMinMover(std::string const & name);
	SaneMinMover(
		core::kinematics::MoveMapOP movemap_in,
		ScoreFunctionOP scorefxn_in,
		MinimizerOptionsOP min_options_in,
		bool cartesian = false
	);

	virtual ~SaneMinMover();
	virtual MoverOP clone() const;

	bool cartesian() const;
	void cartesian( bool setting );

	MinimizerOptionsOP min_options();
	void movemap( core::kinematics::MoveMapOP movemap_in );
	core::kinematics::MoveMapOP movemap();

	void score_function( core::scoring::ScoreFunctionOP scorefxn_in );
	core::scoring::ScoreFunctionOP score_function();

	/// @brief Minimizes the DOFs of pose specified in the MoveMap
	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

private:
	// set reasonable defaults for scorefxn_, movemap_ and min_options_
	void set_defaults_();

	bool cartesian_;
	core::kinematics::MoveMapOP movemap_;
	ScoreFunctionOP scorefxn_;
	MinimizerOptionsOP min_options_;
};

} // moves
} // protocols

#endif

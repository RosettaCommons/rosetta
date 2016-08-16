// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CyclizationMover.hh
/// @brief Header file for CyclizationMover
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_simple_moves_CyclizationMover_hh
#define INCLUDED_protocols_simple_moves_CyclizationMover_hh

// unit headers
#include <protocols/simple_moves/CyclizationMover.fwd.hh>

// protocols headers
#include <protocols/moves/Mover.hh>

// core headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/mm/MMTorsionLibrary.hh>

#include <core/kinematics/MoveMap.hh>

namespace protocols {
namespace simple_moves {

/// @details This mover performas a variaty of functions related to cyclizing a chain in a pose.
class CyclizationMover : public protocols::moves::Mover {
public:

	// ctor
	CyclizationMover( core::Size chain_to_cyclize, bool add_constraints, bool minimize, core::Size minimization_rebuild_rounds );

	// ctor
	CyclizationMover( core::Size chain_to_cyclize, bool add_constraints, bool minimize, core::Size minimization_rebuild_rounds,
		core::scoring::ScoreFunctionOP score_fxn, core::kinematics::MoveMapOP move_map );

	// dtor
	virtual ~CyclizationMover(){}

	// mover interface
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "CyclizationMover"; }

private:
	// mover specific
	void setup_connections( core::pose::Pose & pose );
	void setup_constraints( core::pose::Pose & pose );
	void setup_scorefunction();
	void setup_minimizer( core::pose::Pose & pose );
	void minimize_rebuild( core::pose::Pose & pose );

private:
	// the chain number to cyclize
	core::Size chain_to_cyclize_;

	// N-terminus and C-terminus residue positions
	core::Size nterm_rsd_num_;
	core::Size cterm_rsd_num_;

	// bools used to decide if we should do certain things, functions should be self-evident from name
	bool add_constraints_;
	bool minimize_;

	// rounds of minimization followed by rebuilding residue dependant atom positions to perform
	core::Size minimization_rebuild_rounds_;

	// we are going to need a scorefunction
	//    need setter
	//    need getter
	//    need set at ctor
	//    need default: double check if constraint weights are set to specific values, orther wise set them
	core::scoring::ScoreFunctionOP score_fxn_;

	// we are going to need a move map
	//    need setter
	//    need getter
	//    need set at ctor
	//    need default: all bb and chi free
	core::kinematics::MoveMapOP move_map_;

	/// MMTorsionLibrary for looking up torsion parameters
	core::scoring::mm::MMTorsionLibrary const & mm_torsion_library_;
};

}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_CyclizationMover_HH

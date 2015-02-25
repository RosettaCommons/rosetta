// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/DisulfideInsertionMover.hh
/// @brief Header file for DisulfideInsertionMover
/// @author Orly Marcu ( orly.marcu@mail.huji.ac.il )
/// @date Jan. 12, 2015

#ifndef INCLUDED_protocols_simple_moves_DisulfideInsertionMover_hh
#define INCLUDED_protocols_simple_moves_DisulfideInsertionMover_hh

// unit headers
#include <protocols/simple_moves/DisulfideInsertionMover.fwd.hh>

// protocols headers
#include <protocols/moves/Mover.hh>

// core headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/kinematics/MoveMap.hh>

namespace protocols {
namespace simple_moves {

enum DisulfideCyclizationViability {
	DCV_NOT_CYCLIZABLE,
	DCV_CYCLIZABLE,
	DCV_ALREADY_CYCLIZED
};

/// @detailed a mover that given a receptor peptide pose mutates the peptides edge residues to cysteins,
///           if needed, and enforces disulfide bonding by constrained minimization of the bond and the interaction
class DisulfideInsertionMover : public protocols::moves::Mover {
public:

	// ctor with all parameters
	DisulfideInsertionMover(core::Size const peptide_chain, core::scoring::ScoreFunctionOP scorefxn, core::kinematics::MoveMapOP mm);

	// ctor with default parameters
	DisulfideInsertionMover(core::Size const peptide_chain);

	// dtor
	virtual ~DisulfideInsertionMover(){}

	// mover interface
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "DisulfideInsertionMover"; }

	/// @brief checks if residues next to a putative derived peptide are near enough in space to be mutated to cysteins that might form a cyclic peptide.
	static DisulfideCyclizationViability determine_cyclization_viability(
		core::pose::Pose const & partner_pose,
		core::Size const n_putative_cyd,
		core::Size const c_putative_cyd);

	void set_scorefxn ( core::scoring::ScoreFunctionOP score_function) { scorefxn_ = score_function; }
	core::scoring::ScoreFunctionOP get_scorefxn () { return scorefxn_;}

	void set_movemap (core::kinematics::MoveMapOP mm) {movemap_ = mm;}
	core::kinematics::MoveMapOP get_movemap () {return movemap_;}

	core::Size get_peptide_chain () {return peptide_chain_num_;}

private:
	void setup_constraints(core::pose::Pose & peptide_receptor_pose,
		core::conformation::Residue lower_cys,
		core::conformation::Residue upper_cys,
		core::Size n_cyd_two_chain_position,
		core::Size c_cyd_two_chain_position);


	core::Size peptide_chain_num_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::kinematics::MoveMapOP movemap_;
};

}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_DisulfideInsertionMover_HH

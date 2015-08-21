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

// utility headers
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace simple_moves {

enum DisulfideCyclizationViability {
	DCV_NOT_CYCLIZABLE,
	DCV_CYCLIZABLE,
	DCV_ALREADY_CYCLIZED
};

/// @details a mover that given a receptor peptide pose mutates the peptides edge residues to cysteins,
///           if needed, and enforces disulfide bonding by constrained minimization of the bond and the interaction
class DisulfideInsertionMover : public protocols::moves::Mover {
public:

	// ctor
	DisulfideInsertionMover(core::Size const peptide_chain,
		core::scoring::ScoreFunctionOP scorefxn = NULL, core::kinematics::MoveMapOP mm = NULL,
		bool const is_cyd_res_at_termini = true,
		core::Size const n_cyd_seqpos = 0, core::Size const c_cyd_seqpos = 0);

	// cctor
	DisulfideInsertionMover( DisulfideInsertionMover const & );

	// dtor
	virtual ~DisulfideInsertionMover(){}

	// clone
	virtual protocols::moves::MoverOP clone() const;

	// mover interface
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "DisulfideInsertionMover"; }

	/// @brief checks if residues next to a putative derived peptide are near enough in space to be mutated to cysteins that might form a cyclic peptide.
	static DisulfideCyclizationViability determine_cyclization_viability(
		core::pose::Pose const & partner_pose,
		core::Size const n_putative_cyd_res,
		core::Size const c_putative_cyd_res);

	// RosettaScripts implementation
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	void set_scorefxn ( core::scoring::ScoreFunctionOP const score_function) { scorefxn_ = score_function; }
	core::scoring::ScoreFunctionOP get_scorefxn() const { return scorefxn_;}

	void set_movemap (core::kinematics::MoveMapOP const mm) {movemap_ = mm;}
	core::kinematics::MoveMapOP get_movemap() const {return movemap_;}

	void set_peptide_chain(core::Size const value) { peptide_chain_num_ = value; }
	core::Size get_peptide_chain () const {return peptide_chain_num_;}

	void set_cyd_seqpos(core::Size const n_cyd_seqpos, core::Size const c_cyd_seqpos) { n_cyd_seqpos_ = n_cyd_seqpos; c_cyd_seqpos_ = c_cyd_seqpos; is_cyd_res_at_termini_ = false; }
	void set_cyd_res_at_termini() { is_cyd_res_at_termini_ = true; }
	bool get_is_cyd_res_at_termini() const { return is_cyd_res_at_termini_; }
	core::Size get_n_cyd_seqpos() const { return n_cyd_seqpos_; }
	core::Size get_c_cyd_seqpos() const { return c_cyd_seqpos_; }

	void set_constraint_weight(core::Real const value) { constraint_weight_ = value; }
	core::Real get_constraint_weight() const { return constraint_weight_; }

private:
	/// @brief adds angle, dihedral angle and atom-pair constraints to the pose
	/// Based on code by Nir London
	static void setup_constraints(core::pose::Pose & peptide_receptor_pose,
		core::conformation::Residue lower_cys,
		core::conformation::Residue upper_cys,
		core::Size n_cyd_two_chain_position,
		core::Size c_cyd_two_chain_position);

	/// @brief the chain number for the chain to introduce the Cys mutations to
	core::Size peptide_chain_num_;

	/// @brief the score function to use in all operations; if constraint_weight > 0
	///        this function may be modified to include constraints.
	core::scoring::ScoreFunctionOP scorefxn_;

	/// @brief the movemap for minimization upon call to rebuild_disulfide() (after mutation)
	core::kinematics::MoveMapOP movemap_;

	/// @brief together with c_cyd_seqpos_, the residue numbers to be mutates into Cys
	///        and modeled as having a disulfide bond.
	core::Size n_cyd_seqpos_;
	core::Size c_cyd_seqpos_;

	/// @brief when > 0, applies angle, dihedral angle and atom-pair constraints
	///        using this weight in the score function, to the newly formed disulfide
	///        bond before repacking and minimization.
	core::Real constraint_weight_;

	/// @brief whether the cys seqpos' should be set to the peptide_chain termini.
	/// If this is true, any setting to c_cyd_seqpos_ and n_cyd_seqpos_ is discarded
	/// upon a call to apply().
	bool is_cyd_res_at_termini_;
};

}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_DisulfideInsertionMover_HH


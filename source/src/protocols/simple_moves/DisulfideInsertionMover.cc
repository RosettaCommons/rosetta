// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/DisulfideInsersionMover.cc
/// @brief a mover that closes a receptor bound peptide by an added disulfide bond
/// @author Orly Marcu ( orlymarcu@mail.huji.ac.il )
/// @date Jan. 12, 2014

// unit headers
#include <protocols/simple_moves/DisulfideInsertionMover.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/util/disulfide_util.hh>


//Utility headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
//#include <basic/options/option.hh>


namespace protocols {
namespace simple_moves {

static thread_local basic::Tracer TR( "protocols.simple_moves.DisulfideInsertionMover" );

/// @brief constructor with all parameters
DisulfideInsertionMover::DisulfideInsertionMover(core::Size const peptide_chain, core::scoring::ScoreFunctionOP scorefxn, core::kinematics::MoveMapOP mm) :
	protocols::moves::Mover("DisulfideInsertionMover"),	peptide_chain_num_(peptide_chain),
	scorefxn_(scorefxn),
	movemap_( mm )
{
}

/// @brief constructor with default parameters
DisulfideInsertionMover::DisulfideInsertionMover(core::Size peptide_chain) :
	protocols::moves::Mover("DisulfideInsertionMover"),
	peptide_chain_num_(peptide_chain),
	scorefxn_(core::scoring::get_score_function()),
	movemap_( new core::kinematics::MoveMap)
{
	movemap_->set_chi(true);
	movemap_->set_jump(false);
}

/// Finds out if the residues next to the peptide chain (assumed second in pose)
/// are in a specified distance range (3-5 Ang, unless one of the residues is Gly,
/// in which case 4.5-6.7 Ang) so that the peptide they wrap can be stabilized by
/// cyclization.
///
/// @param partner_pose the protein from which peptides are derived
/// @param n_putative_cyd position in the partner pose where a putative cysteine immediatly before a derived peptide might exist
/// @param c_putative_cyd position in the partner pose where a putative cysteine immediatly after a derived peptide might exist
DisulfideCyclizationViability
DisulfideInsertionMover::determine_cyclization_viability(
		core::pose::Pose const & partner_pose,
		core::Size const n_putative_cyd, core::Size const c_putative_cyd) {

	//lower and upper boundries for atom distances that indicate closability by a disulfide
	core::Length min_distance;
	core::Length max_distance;

	std::string atomToCheck;
	if ((partner_pose.residue(n_putative_cyd).name1() == 'G' ) ||( partner_pose.residue( c_putative_cyd ).name1() == 'G' )) {
		atomToCheck	= "CA";
		min_distance = 4.5;
		max_distance = 6.5;
	}
	else {
		atomToCheck = "CB";
		min_distance = 3;
		max_distance = 5;
	}
	core::Length distance =	partner_pose.residue( n_putative_cyd ).xyz( atomToCheck ).distance(partner_pose.residue( c_putative_cyd ).xyz( atomToCheck ));
	bool already_bonded_n_cys (partner_pose.conformation().residue(n_putative_cyd).has_variant_type( core::chemical::DISULFIDE ));
	bool already_bonded_c_cys (partner_pose.conformation().residue(c_putative_cyd).has_variant_type( core::chemical::DISULFIDE ));
	// both don't form disulfide bonds and so their potential can be estimated according to distance alone
	if (!already_bonded_n_cys && !already_bonded_c_cys) {
		if (distance <= max_distance && distance >= min_distance) {
			return DCV_CYCLIZABLE;
		}
	}
	else if (already_bonded_n_cys && already_bonded_c_cys) {
		// both cysteins are bonded but not neccessarily to each other
		core::Size lower_SG_id = partner_pose.residue(n_putative_cyd).atom_index( "SG" );
		core::Size cys_SG_connect_id =  partner_pose.residue(n_putative_cyd).type().residue_connection_id_for_atom( lower_SG_id );
		if (partner_pose.residue(n_putative_cyd).residue_connection_partner(cys_SG_connect_id) == c_putative_cyd) {
			return DCV_ALREADY_CYCLIZED;
		}
	}
	return DCV_NOT_CYCLIZABLE;
}

void
DisulfideInsertionMover::apply( core::pose::Pose & peptide_receptor_pose )
{

	core::Size n_cyd_two_chain_position = peptide_receptor_pose.conformation().chain_begin(peptide_chain_num_);
	core::Size c_cyd_two_chain_position = peptide_receptor_pose.total_residue();

	// eliminate cases where a disulfide should not be formed since the residues already form a disulfide (closability==2)
	// in that case we will only want to optimize it using the rebuild_disulfide function
	if (determine_cyclization_viability(peptide_receptor_pose, n_cyd_two_chain_position, c_cyd_two_chain_position) == DCV_CYCLIZABLE) {
		core::conformation::form_disulfide(peptide_receptor_pose.conformation(), n_cyd_two_chain_position, c_cyd_two_chain_position);
	}

	core::pose::Pose original_orientation(peptide_receptor_pose);
	core::conformation::Residue lower_cys = peptide_receptor_pose.residue( n_cyd_two_chain_position);
	core::conformation::Residue upper_cys = peptide_receptor_pose.residue( c_cyd_two_chain_position);
	setup_constraints(peptide_receptor_pose, lower_cys, upper_cys, n_cyd_two_chain_position, c_cyd_two_chain_position);
	// TODO : No matter if movemap is provided or initialized, peptide backbone is free to move
	// probably a must for the disulfide to be properly rebuilt, but need to think about this forced DOF
	movemap_->set_bb_true_range(n_cyd_two_chain_position, c_cyd_two_chain_position);
	core::util::rebuild_disulfide(peptide_receptor_pose, n_cyd_two_chain_position, c_cyd_two_chain_position, NULL, scorefxn_, movemap_, scorefxn_);
	scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0 );
	scorefxn_->set_weight( core::scoring::dihedral_constraint, 0 );
	scorefxn_->set_weight( core::scoring::angle_constraint, 0 );

}


/// @brief Setup the score function with the appropriate weights on the constraints.
/// If user has provided a score function the and set the weights to a specific value,
/// that value is kept. If weight is not set, prints a warning.
void
DisulfideInsertionMover::setup_constraints( core::pose::Pose & peptide_receptor_pose,
											core::conformation::Residue lower_cys,
											core::conformation::Residue upper_cys,
											core::Size n_cyd_two_chain_position,
											core::Size c_cyd_two_chain_position)
{
	core::scoring::constraints::ConstraintSetOP cys_cst_set( new core::scoring::constraints::ConstraintSet() );
	core::id::AtomID n_ca_atom_id( lower_cys.atom_index("CA"),n_cyd_two_chain_position  );
	core::id::AtomID n_cb_atom_id( lower_cys.atom_index("CB"), n_cyd_two_chain_position );
	core::id::AtomID n_sg_atom_id( lower_cys.atom_index("SG"), n_cyd_two_chain_position );
	core::id::AtomID c_ca_atom_id( upper_cys.atom_index("CA"), c_cyd_two_chain_position );
	core::id::AtomID c_cb_atom_id( upper_cys.atom_index("CB"), c_cyd_two_chain_position );
	core::id::AtomID c_sg_atom_id( upper_cys.atom_index("SG"), c_cyd_two_chain_position );

	cys_cst_set->add_constraint(
				core::scoring::constraints::AtomPairConstraintCOP(new core::scoring::constraints::AtomPairConstraint(
																  n_sg_atom_id,
																  c_sg_atom_id,
																  core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 2.02 /*anchor*/, 0.35 /*stdev*/)))));
	cys_cst_set->add_constraint(
				core::scoring::constraints::DihedralConstraintCOP(new core::scoring::constraints::DihedralConstraint(
																  n_ca_atom_id,
																  n_cb_atom_id,
																  n_sg_atom_id,
																  c_sg_atom_id,
																  core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( -0.558 /*anchor*/, 0.1 /*stdev*/) ))));
	cys_cst_set->add_constraint(
				core::scoring::constraints::DihedralConstraintCOP(new core::scoring::constraints::DihedralConstraint(
																  c_ca_atom_id,
																  c_cb_atom_id,
																  c_sg_atom_id,
																  n_sg_atom_id,
																  core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( -2 /*anchor*/, 0.1 /*stdev*/) ))));
	cys_cst_set->add_constraint(
			     core::scoring::constraints::AngleConstraintCOP(new core::scoring::constraints::AngleConstraint(
		    		 	 	 	 	 	 	 	 	 	 		n_cb_atom_id,
   																n_sg_atom_id,
																c_sg_atom_id,
																core::scoring::func::HarmonicFuncOP(new core::scoring::func::HarmonicFunc( 1.804 /*anchor*/, 0.08 /*stdev*/) )) ));
	cys_cst_set->add_constraint(
				core::scoring::constraints::AngleConstraintCOP(new core::scoring::constraints::AngleConstraint(
															   c_cb_atom_id,
															   c_sg_atom_id,
															   n_sg_atom_id,
															   core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.804 /*anchor*/, 0.08 /*stdev*/) )) ));
	peptide_receptor_pose.constraint_set( cys_cst_set );
	scorefxn_->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	scorefxn_->set_weight( core::scoring::dihedral_constraint, 1.0 );
	scorefxn_->set_weight( core::scoring::angle_constraint, 1.0 );

}


}//simple_moves
}//protocols

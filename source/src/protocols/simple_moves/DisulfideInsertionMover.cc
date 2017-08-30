// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/DisulfideInsertionMover.cc
/// @brief a mover that closes a receptor bound peptide by an added disulfide bond
/// @author Orly Marcu ( orly.marcu@mail.huji.ac.il )
/// @author Yuval Sedan
/// @date Jan. 12, 2014

// unit headers
#include <protocols/simple_moves/DisulfideInsertionMover.hh>

// RosettaScripts headers
#include <protocols/rosetta_scripts/util.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/util/disulfide_util.hh>
#include <core/io/pdb/build_pose_as_is.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/DisulfideInsertion.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

//remove after removing prints
#include <utility/io/ozstream.hh>
#include <core/io/pdb/pdb_writer.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/simple_moves/MinMover.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/simple_moves/DisulfideInsertionMoverCreator.hh>

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.DisulfideInsertionMover" );



/// @brief empty constructor - member variables initialized from the options system
DisulfideInsertionMover::DisulfideInsertionMover() :
	protocols::moves::Mover("DisulfideInsertionMover")
{
	peptide_chain_num_= basic::options::option[ basic::options::OptionKeys::DisulfideInsertion::peptide_chain ]();
	n_cyd_seqpos_ = basic::options::option[ basic::options::OptionKeys::DisulfideInsertion::n_cyd_seqpos ]();
	c_cyd_seqpos_ = basic::options::option[ basic::options::OptionKeys::DisulfideInsertion::n_cyd_seqpos ]();
	max_dslf_pot_ = basic::options::option[ basic::options::OptionKeys::DisulfideInsertion::max_dslf_pot ]();
	max_dslf_energy_ = basic::options::option [ basic::options::OptionKeys::DisulfideInsertion::max_dslf_energy ]();
	min_dist_multiplier_ = basic::options::option[ basic::options::OptionKeys::DisulfideInsertion::min_dslf_dist_multiplier ]();
	max_dist_multiplier_ = basic::options::option[ basic::options::OptionKeys::DisulfideInsertion::max_dslf_dist_multiplier ]();
	constraint_weight_ = basic::options::option[ basic::options::OptionKeys::DisulfideInsertion::constraint_weight ]();

	if ( movemap_ == nullptr ) {
		movemap_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	}

	if ( scorefxn_ == nullptr ) {
		scorefxn_ = core::scoring::get_score_function();
	}
}

DisulfideInsertionMover::DisulfideInsertionMover(DisulfideInsertionMover const &other) :
	protocols::moves::Mover("DisulfideInsertionMover") {
	set_peptide_chain( other.get_peptide_chain() );
	set_n_cyd_seqpos(other.get_n_cyd_seqpos() );
	set_c_cyd_seqpos(other.get_c_cyd_seqpos());
	// TODO : should we clone this? simply calling ->clone() didn't work for some reason (at least, get_name() on the new scorefxn returned an empty string)
	set_scorefxn(other.get_scorefxn());
	set_movemap(other.get_movemap());
	set_constraint_weight(other.get_constraint_weight());
}

// @brief clone function that returns a deep copy of the mover
protocols::moves::MoverOP DisulfideInsertionMover::clone() const {
	return DisulfideInsertionMoverOP(new DisulfideInsertionMover(*this));
}



/// Finds out if the residues next to the peptide chain (assumed second in pose)
/// are in a specified distance range (3-5 Ang, unless one of the residues is Gly,
/// in which case 4.5-6.5 Ang) so that the peptide they wrap can be stabilized by
/// cyclization.
///
/// @param pose the protein from which peptides are derived
/// @param n_putative_cyd position in the partner pose where a putative cysteine immediatly before a derived peptide might exist
/// @param c_putative_cyd position in the partner pose where a putative cysteine immediatly after a derived peptide might exist

DisulfideCyclizationViability
DisulfideInsertionMover::determine_cyclization_viability(
	core::pose::Pose const & pose,
	core::Size const n_putative_cyd_idx, core::Size const c_putative_cyd_idx) {

	if ( n_putative_cyd_idx < pose.conformation().chain_begin(1) || c_putative_cyd_idx > pose.total_residue() ) {
		TR << "one or more of the requested positions " << n_putative_cyd_idx << " and " << c_putative_cyd_idx << "do not exist in this pose" << std::endl;
		return DCV_NOT_CYCLIZABLE;
	}
	// assume it is safe to hold a const reference to the Residues, since we are given
	// a const Pose
	core::conformation::Residue const & n_putative_cyd( pose.residue(n_putative_cyd_idx) );
	core::conformation::Residue const & c_putative_cyd( pose.residue(c_putative_cyd_idx) );

	if ( n_putative_cyd.chain()!=c_putative_cyd.chain() ) {
		TR << "Requested positions not found in this same chain. Cannot form intra-chain disulfide bond" << std::endl;
		return DCV_NOT_CYCLIZABLE;
	}

	if ( !n_putative_cyd.is_protein() || !c_putative_cyd.is_protein() ) {
		return DCV_NOT_CYCLIZABLE;
	}

	//lower and upper boundries for atom distances that indicate closability by a disulfide
	core::Length min_distance;
	core::Length max_distance;

	// first, check if the requested residues are already disulfide-bonded
	// (implicitly, if the residue has a DISULFIDE variant type, we assume it is CYS)
	bool already_bonded_n_cys (n_putative_cyd.has_variant_type( core::chemical::DISULFIDE ));
	bool already_bonded_c_cys (c_putative_cyd.has_variant_type( core::chemical::DISULFIDE ));

	if ( !already_bonded_n_cys && !already_bonded_c_cys ) {
		// both residues don't already have a disulfide bond

		// estimate potential for cyclization by distance (CB-CB)
		// glycines don't have CBs, so we use CA-CA distance
		std::string atom_to_check;
		if ( (n_putative_cyd.name1() == 'G') ||
				(c_putative_cyd.name1() == 'G') ) {
			atom_to_check = "CA";
			min_distance = 4.5*min_dist_multiplier_;
			max_distance = 6.5*max_dist_multiplier_;
		} else {
			atom_to_check = "CB";
			min_distance = 3*min_dist_multiplier_;
			max_distance = 5*max_dist_multiplier_;
		}

		core::Length distance = n_putative_cyd.xyz( atom_to_check ).distance(c_putative_cyd.xyz( atom_to_check ));

		if ( distance > max_distance || distance < min_distance ) {
			TR.Debug << " requested positions " << n_putative_cyd_idx << " and " << c_putative_cyd_idx << " are not in the required distance range (distance = " << distance << ", max = " << max_distance << ", min = " << min_distance << " )"  << std::endl;
			return DCV_NOT_CYCLIZABLE;
		} else {
			TR.Debug << " requested positions " << n_putative_cyd_idx << " and " << c_putative_cyd_idx << " ARE in the required distance range (distance = " << distance << ", max = " << max_distance << ", min = " << min_distance << " )"  << std::endl;
		}

		core::scoring::disulfides::DisulfideMatchingPotential disulfPot;
		core::Energy match_t = 0.0;
		core::Energy match_r = 0.0;
		core::Energy match_rt = 0.0;
		// TODO check this for dcys too?
		disulfPot.score_disulfide( n_putative_cyd, c_putative_cyd, match_t, match_r, match_rt, false /*mirror score for D-cys potential*/ );
		TR << " positions " << n_putative_cyd_idx << " and " << c_putative_cyd_idx << " have a match potential of " << match_rt << " (cutoff = " << max_dslf_pot_ << " )" << std::endl;
		if ( match_rt <= max_dslf_pot_ ) {
			return DCV_CYCLIZABLE;
		}
		// the residues are not at the right distance from each other
		TR.Debug << " requested positions " << n_putative_cyd_idx << " and " << c_putative_cyd_idx << " don't have the required match potential " << std::endl;
		return DCV_NOT_CYCLIZABLE;
	} else if ( already_bonded_n_cys ) {
		// the N-terminal is a disulfide-bonded cysteine.
		// check what the N-terminal cysteine is bonded to
		core::Size n_cys_SG_id = n_putative_cyd.atom_index( "SG" );
		core::Size n_cys_SG_connect_id =  n_putative_cyd.type().residue_connection_id_for_atom( n_cys_SG_id );
		core::Size resi_bound_to_n_cyd = n_putative_cyd.residue_connection_partner(n_cys_SG_connect_id);

		// if the bond is between the N- and C- residues, meaning that the peptide is already cyclic
		if ( resi_bound_to_n_cyd == c_putative_cyd_idx ) {
			return DCV_ALREADY_CYCLIZED;
		} else {
			// the N-terminal cysteine is bonded to another residue
			// in this scenario, we don't replace it, since it is
			// less likely it will form a disulfide bond with a cysteine we
			// will introduce.
			TR.Debug << " One or more of the requested positions " << n_putative_cyd_idx << " and " << c_putative_cyd_idx << " already involved in a disulfide bond" << std::endl;
			return DCV_NOT_CYCLIZABLE;
		}

	}

	// the C-terminal cysteine has a disulfide bond, but it is not to the N-terminal cysteine
	// see above for our considerations.
	TR.Debug << " One or more of the requested positions " << n_putative_cyd_idx << " and " << c_putative_cyd_idx << " already involved in a disulfide bond" << std::endl;
	return DCV_NOT_CYCLIZABLE;
}

void
DisulfideInsertionMover::apply( core::pose::Pose & peptide_receptor_pose )
{
	set_last_move_status(protocols::moves::FAIL_RETRY);
	core::Size this_pose_n_cyd;
	core::Size this_pose_c_cyd;

	if ( get_n_cyd_seqpos()==0 ) {
		core::Size peptide_start_pos = peptide_receptor_pose.conformation().chain_begin(peptide_chain_num_);
		this_pose_n_cyd = peptide_start_pos;
	} else {
		this_pose_n_cyd = n_cyd_seqpos_;
	}
	if ( get_c_cyd_seqpos()==0 ) {
		core::Size peptide_end_pos =  peptide_receptor_pose.conformation().chain_end(peptide_chain_num_);
		this_pose_c_cyd = peptide_end_pos;
	} else {
		this_pose_c_cyd = c_cyd_seqpos_;
	}

	protocols::simple_moves::DisulfideCyclizationViability cyclizable = determine_cyclization_viability(peptide_receptor_pose,this_pose_n_cyd, this_pose_c_cyd);

	if ( cyclizable == DCV_NOT_CYCLIZABLE ) {
		return;
	}

	peptide_receptor_pose.conformation().detect_disulfides();
	core::scoring::EnergyMap original_pose_emap (peptide_receptor_pose.energies().total_energies());
	core::Real baseline_dslf_energy = original_pose_emap.get(core::scoring::score_type_from_name("dslf_fa13"));

	// eliminate cases where a disulfide should not be formed since the residues already form a disulfide (DCV_ALREADY_CYCLIZED)
	// in that case we will only want to optimize it using the rebuild_disulfide function
	if ( cyclizable == DCV_CYCLIZABLE ) {
		core::conformation::form_disulfide(peptide_receptor_pose.conformation(), this_pose_n_cyd, this_pose_c_cyd);
	}

	core::conformation::Residue lower_cys = peptide_receptor_pose.residue( this_pose_n_cyd);
	core::conformation::Residue upper_cys = peptide_receptor_pose.residue( this_pose_c_cyd);

	if ( constraint_weight_ > 0 ) { // save some code-running if no constraint it set
		scorefxn_->set_weight( core::scoring::atom_pair_constraint, constraint_weight_ );
		scorefxn_->set_weight( core::scoring::dihedral_constraint, constraint_weight_ );
		scorefxn_->set_weight( core::scoring::angle_constraint, constraint_weight_ );

		DisulfideInsertionMover::setup_constraints(peptide_receptor_pose, lower_cys, upper_cys, this_pose_n_cyd, this_pose_c_cyd);
	}

	// TODO : No matter if movemap is provided or initialized, peptide backbone is free to move
	// probably a must for the disulfide to be properly rebuilt, but need to think about this forced DOF
	movemap_->set_bb_true_range(this_pose_n_cyd, this_pose_c_cyd);
	core::util::rebuild_disulfide(peptide_receptor_pose, this_pose_n_cyd, this_pose_c_cyd,
		/*packer_task=*/nullptr,
		/*packer_score=*/scorefxn_,
		/*mm=*/movemap_,
		/*minimizer_score=*/scorefxn_);

	// add our dfpmin minimization - increases number of cases where successful closure is achieved compared to lbfgs_armijo_nonmonotone
	protocols::simple_moves::MinMover minimizer( movemap_, scorefxn_, "dfpmin_armijo_atol", 0.01 /*tolerance*/, true /*nb_list*/ );
	minimizer.apply( peptide_receptor_pose );

	if ( constraint_weight_ > 0 ) { // save some code-running if no constraint it set
		scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0 );
		scorefxn_->set_weight( core::scoring::dihedral_constraint, 0 );
		scorefxn_->set_weight( core::scoring::angle_constraint, 0 );
	}


	core::scoring::EnergyMap cyclic_pose_emap (peptide_receptor_pose.energies().total_energies());
	core::Real const cyc_pose_dslf_energy = cyclic_pose_emap.get(core::scoring::score_type_from_name("dslf_fa13"));
	core::Real const dis_energy_change = cyc_pose_dslf_energy - baseline_dslf_energy;

	TR << "original peptide disulfide energy: " << baseline_dslf_energy << " cyclic peptide disulfide energy: " << cyc_pose_dslf_energy << std::endl;
	TR << "change in disulfide energy upon peptide disulfide bridge formation: " << dis_energy_change << std::endl;
	if ( dis_energy_change < max_dslf_energy_ ) {
		set_last_move_status(protocols::moves::MS_SUCCESS);
	}
}



/// @brief Setup the score function with the appropriate weights on the constraints.
///        If user has provided a score function and set the weights to a specific value,
///        that value is kept. If weight is not set, prints a warning.
///        This is based on a function written by Nir London.
void
DisulfideInsertionMover::setup_constraints( core::pose::Pose & peptide_receptor_pose,
	core::conformation::Residue lower_cys,
	core::conformation::Residue upper_cys,
	core::Size n_cyd_seqpos,
	core::Size c_cyd_seqpos)
{
	core::scoring::constraints::ConstraintSetOP cys_cst_set( new core::scoring::constraints::ConstraintSet() );
	core::id::AtomID n_ca_atom_id( lower_cys.atom_index("CA"),n_cyd_seqpos  );
	core::id::AtomID n_cb_atom_id( lower_cys.atom_index("CB"), n_cyd_seqpos );
	core::id::AtomID n_sg_atom_id( lower_cys.atom_index("SG"), n_cyd_seqpos );
	core::id::AtomID c_ca_atom_id( upper_cys.atom_index("CA"), c_cyd_seqpos );
	core::id::AtomID c_cb_atom_id( upper_cys.atom_index("CB"), c_cyd_seqpos );
	core::id::AtomID c_sg_atom_id( upper_cys.atom_index("SG"), c_cyd_seqpos );

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

}


// RosettaScripts implementation
void DisulfideInsertionMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose ) {

	if ( tag->hasOption("scorefxn") ) {
		set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data, "" ) );
	}

	// TODO : check if such a check is neccesary: if ( tag->hasOption("MoveMap") ) {
	protocols::rosetta_scripts::parse_movemap( tag, pose, movemap_ );

	if ( tag->hasOption("n_cyd") ) {
		set_n_cyd_seqpos(tag->getOption<core::Size>("n_cyd"));
	}

	if ( tag->hasOption("c_cyd") ) {
		set_c_cyd_seqpos(tag->getOption<core::Size>("c_cyd"));
	}

	if ( tag->hasOption("chain") ) {
		set_peptide_chain(tag->getOption<core::Size>("chain"));
	}

	if ( tag->hasOption("constraint_weight") ) {
		set_constraint_weight(tag->getOption<core::Real>("constraint_weight"));
	}

	if ( tag->hasOption("max_dslf_pot") ) {
		set_max_dslf_pot(tag->getOption<core::Real>("max_dslf_pot"));
	}

}

std::string DisulfideInsertionMover::get_name() const {
	return mover_name();
}

std::string DisulfideInsertionMover::mover_name() {
	return "DisulfideInsertion";
}

void DisulfideInsertionMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function(attlist);

	XMLSchemaSimpleSubelementList ssl;
	rosetta_scripts::append_subelement_for_parse_movemap_w_datamap(xsd, ssl);

	attlist + XMLSchemaAttribute(
		"n_cyd", xsct_non_negative_integer,
		"XSD XRW: TO DO specifiy n_cyd and c_cyd or neither");

	attlist + XMLSchemaAttribute(
		"c_cyd", xsct_non_negative_integer,
		"XSD XRW: TO DO specifiy n_cyd and c_cyd or neither");

	attlist + XMLSchemaAttribute(
		"chain", xsct_non_negative_integer,
		"XSD XRW: TO DO");

	attlist + XMLSchemaAttribute(
		"constraint_weight", xsct_real,
		"XSD XRW: TO DO");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"XSD XRW: TO DO",
		attlist, ssl);
}

std::string DisulfideInsertionMoverCreator::keyname() const {
	return DisulfideInsertionMover::mover_name();
}

protocols::moves::MoverOP
DisulfideInsertionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DisulfideInsertionMover );
}

void DisulfideInsertionMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DisulfideInsertionMover::provide_xml_schema( xsd );
}


}//simple_moves
}//protocols

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
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/simple_moves/DisulfideInsertionMoverCreator.hh>

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.DisulfideInsertionMover" );

/// @brief constructor with all parameters
DisulfideInsertionMover::DisulfideInsertionMover(core::Size const peptide_chain,
	core::scoring::ScoreFunctionOP scorefxn, core::kinematics::MoveMapOP mm,
	bool const is_cyd_res_at_termini,
	core::Size const n_cyd_seqpos, core::Size const c_cyd_seqpos) :
	protocols::moves::Mover("DisulfideInsertionMover"),
	peptide_chain_num_( peptide_chain ),
	scorefxn_(std::move( scorefxn )),
	movemap_(std::move( mm )),
	n_cyd_seqpos_(n_cyd_seqpos),
	c_cyd_seqpos_(c_cyd_seqpos),
	is_cyd_res_at_termini_(is_cyd_res_at_termini)
{
	if ( movemap_ == nullptr ) {
		movemap_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	}
	if ( scorefxn_ == nullptr ) {
		scorefxn_ = core::scoring::get_score_function();
	}
	constraint_weight_ = basic::options::option[ basic::options::OptionKeys::run::insert_disulfide_constraint_weight ]();
}

DisulfideInsertionMover::DisulfideInsertionMover(DisulfideInsertionMover const &other) :
	protocols::moves::Mover("DisulfideInsertionMover") {
	set_peptide_chain( other.get_peptide_chain() );
	if ( other.get_is_cyd_res_at_termini() ) {
		set_cyd_res_at_termini();
	} else {
		set_cyd_seqpos( other.get_n_cyd_seqpos(), other.get_c_cyd_seqpos() );
	}
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
	core::Size const n_putative_cyd_index, core::Size const c_putative_cyd_index) {

	// assume it is safe to hold a const reference to the Residues, since we are given
	// a const Pose
	core::conformation::Residue const & n_putative_cyd( pose.residue(n_putative_cyd_index) );
	core::conformation::Residue const & c_putative_cyd( pose.residue(c_putative_cyd_index) );

	//lower and upper boundries for atom distances that indicate closability by a disulfide
	core::Length min_distance;
	core::Length max_distance;

	if ( !n_putative_cyd.is_protein() || !c_putative_cyd.is_protein() ) {
		return DCV_NOT_CYCLIZABLE;
	}

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
			min_distance = 4.5;
			max_distance = 6.5;
		} else {
			atom_to_check = "CB";
			min_distance = 3;
			max_distance = 5;
		}

		core::Length distance = n_putative_cyd.xyz( atom_to_check ).distance(c_putative_cyd.xyz( atom_to_check ));

		if ( distance <= max_distance && distance >= min_distance ) {
			return DCV_CYCLIZABLE;
		}

		// the residues are not at the right distance from each other
		return DCV_NOT_CYCLIZABLE;
	} else if ( already_bonded_n_cys ) {
		// the N-terminal is a disulfide-bonded cysteine.

		// check what the N-terminal cysteine is bonded to
		core::Size n_cys_SG_id = n_putative_cyd.atom_index( "SG" );
		core::Size n_cys_SG_connect_id =  n_putative_cyd.type().residue_connection_id_for_atom( n_cys_SG_id );
		core::Size resi_bound_to_n_cyd = n_putative_cyd.residue_connection_partner(n_cys_SG_connect_id);

		// if the bond is between the N- and C- residues, meaning that the peptide is already cyclic
		if ( resi_bound_to_n_cyd == c_putative_cyd_index ) {
			return DCV_ALREADY_CYCLIZED;
		} else {
			// the N-terminal cysteine is bonded to another residue
			// in this scenario, we don't replace it, since it is
			// less likely it will form a disulfide bond with a cysteine we
			// will introduce.
			return DCV_NOT_CYCLIZABLE;
		}

	}

	// the C-terminal cysteine has a disulfide bond, but it is not to the N-terminal cysteine
	// see above for our considerations.
	return DCV_NOT_CYCLIZABLE;
}

void
DisulfideInsertionMover::apply( core::pose::Pose & peptide_receptor_pose )
{
	set_last_move_status(protocols::moves::FAIL_RETRY);
	core::Size peptide_start_pos = peptide_receptor_pose.conformation().chain_begin(peptide_chain_num_);
	core::Size peptide_end_pos =  peptide_receptor_pose.conformation().chain_end(peptide_chain_num_);

	if ( is_cyd_res_at_termini_ ) {
		n_cyd_seqpos_ = peptide_start_pos;
		c_cyd_seqpos_ = peptide_end_pos;
	}
	protocols::simple_moves::DisulfideCyclizationViability cyclizable = determine_cyclization_viability(peptide_receptor_pose, n_cyd_seqpos_, c_cyd_seqpos_);

	if ( cyclizable == DCV_NOT_CYCLIZABLE ) {
		return;
	}

	// a copy of the pose to compare energies to after disulfide insertion step
	//core::pose::PoseOP original_peptide (new core::pose::Pose (peptide_receptor_pose, n_cyd_seqpos_, c_cyd_seqpos_) );

	// eliminate cases where a disulfide should not be formed since the residues already form a disulfide (DCV_ALREADY_CYCLIZED)
	// in that case we will only want to optimize it using the rebuild_disulfide function
	if ( cyclizable == DCV_CYCLIZABLE ) {
		core::conformation::form_disulfide(peptide_receptor_pose.conformation(), n_cyd_seqpos_, c_cyd_seqpos_);
	}

	core::conformation::Residue lower_cys = peptide_receptor_pose.residue( n_cyd_seqpos_);
	core::conformation::Residue upper_cys = peptide_receptor_pose.residue( c_cyd_seqpos_);

	if ( constraint_weight_ > 0 ) { // save some code-running if no constraint it set
		scorefxn_->set_weight( core::scoring::atom_pair_constraint, constraint_weight_ );
		scorefxn_->set_weight( core::scoring::dihedral_constraint, constraint_weight_ );
		scorefxn_->set_weight( core::scoring::angle_constraint, constraint_weight_ );

		DisulfideInsertionMover::setup_constraints(peptide_receptor_pose, lower_cys, upper_cys, n_cyd_seqpos_, c_cyd_seqpos_);
	}

	// TODO : No matter if movemap is provided or initialized, peptide backbone is free to move
	// probably a must for the disulfide to be properly rebuilt, but need to think about this forced DOF
	movemap_->set_bb_true_range(n_cyd_seqpos_, c_cyd_seqpos_);
	core::util::rebuild_disulfide(peptide_receptor_pose, n_cyd_seqpos_, c_cyd_seqpos_,
		/*packer_task=*/nullptr,
		/*packer_score=*/scorefxn_,
		/*mm=*/movemap_,
		/*minimizer_score=*/scorefxn_);

	if ( constraint_weight_ > 0 ) { // save some code-running if no constraint it set
		scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0 );
		scorefxn_->set_weight( core::scoring::dihedral_constraint, 0 );
		scorefxn_->set_weight( core::scoring::angle_constraint, 0 );
	}

	// evaluate energy of newly formed peptide
	//core::pose::PoseOP cyclic_peptide (new core::pose::Pose (peptide_receptor_pose, n_cyd_seqpos_, c_cyd_seqpos_)  );
	//core::Real cyclic_peptide_energy = (*scorefxn_)(*cyclic_peptide);
	//core::Real original_peptide_energy = (*scorefxn_)(*original_peptide);

	// exclude energetic gain by disulfide when comparing all other aspects of peptide change
	//core::scoring::EnergyMap cyclic_pep_emap (cyclic_peptide->energies().total_energies());
	//core::scoring::EnergyMap original_pep_emap (original_peptide->energies().total_energies());


	//core::Real baseline_dslf_energy = original_pep_emap.get(core::scoring::score_type_from_name("dslf_fa13"));
	//core::Real cyc_pep_dslf_energy = cyclic_pep_emap.get(core::scoring::score_type_from_name("dslf_fa13"));

	//const core::Real REASONABLE_ENERGY_CHANGE = 1.55;
	//core::Real dis_energy_change = cyc_pep_dslf_energy - baseline_dslf_energy;
	//core::Real other_energetic_change = (cyclic_peptide_energy - original_peptide_energy) - dis_energy_change;

	// if both the energy of the new disulfide bond AND of all other terms is not much worse call this a success
	//if (dis_energy_change < REASONABLE_ENERGY_CHANGE && other_energetic_change < REASONABLE_ENERGY_CHANGE) {
	set_last_move_status(protocols::moves::MS_SUCCESS);
	//}
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

	// There are two ways to specify residues to mutate.
	// There are three scenarios for specified options.
	// a. n/c_cyd both not specified => at_termini=true
	// b. only one of n/c_cyd is specified => error
	// c. n/c_cyd specified => use specified n/c_cyd

	bool const n_cyd_specified(tag->hasOption("n_cyd"));
	bool const c_cyd_specified(tag->hasOption("c_cyd"));
	if ( n_cyd_specified && c_cyd_specified ) {
		set_cyd_seqpos(tag->getOption<core::Size>("n_cyd"), tag->getOption<core::Size>("c_cyd"));
	} else if ( n_cyd_specified || c_cyd_specified ) {
		throw utility::excn::EXCN_RosettaScriptsOption("n_cyd and c_cyd must either both be specified or both not specified");
	} else {
		set_cyd_res_at_termini();
	}

	if ( tag->hasOption("chain") ) {
		set_peptide_chain(tag->getOption<core::Size>("chain"));
	}

	if ( tag->hasOption("constraint_weight") ) {
		set_constraint_weight(tag->getOption<core::Real>("constraint_weight"));
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

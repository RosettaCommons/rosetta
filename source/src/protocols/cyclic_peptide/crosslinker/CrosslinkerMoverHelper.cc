// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.cc
/// @brief A base class for helper objects that the CrosslinkerMover uses to set up specific types
/// of threefold linkers.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/Energies.hh>

// Protocols headers
#include <protocols/simple_moves/ClearConstraintsMover.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.CrosslinkerMoverHelper" );

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CrosslinkerMoverHelper::CrosslinkerMoverHelper() //:
//TODO initialize data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
CrosslinkerMoverHelper::CrosslinkerMoverHelper( CrosslinkerMoverHelper const & /*src*/ ) //:
//TODO copy data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
CrosslinkerMoverHelper::~CrosslinkerMoverHelper(){}

//////////////////////
/// Public Methods ///
//////////////////////

/// @brief Given a ResidueSubset with exactly three residues selected, pull out the three indices.
/// @details Overwrites res1, res2, and res3.
void
CrosslinkerMoverHelper::get_sidechain_indices(
	core::select::residue_selector::ResidueSubset const & selection,
	core::Size &res1,
	core::Size &res2,
	core::Size &res3
) const {
	core::Size const nres( selection.size() );
	runtime_assert_string_msg( nres>=3, "Error in protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelper::get_sidechain_indices(): Fewer than three residues are in the pose." );
	res1=0; res2=0; res3=0;
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( selection[i] ) {
			if ( res1==0 ) { res1 = i; }
			else if ( res2==0 ) { res2 = i; }
			else if ( res3==0 ) { res3 = i; }
			else {
				utility_exit_with_message( "Error in protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelper::get_sidechain_indices(): More than three residues were selected." );
			}
		}
	}

	runtime_assert_string_msg( res1>0 && res2>0 && res3>0, "Error in protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelper::get_sidechain_indices(): Fewer than three residues were selected." );
}

/// @brief Determine whether a selection is symmetric.
/// @details Returns true if and only if (a) the pose is symmetric, (b) there are three symmetry copies, and (c) the selected residues are equivalent residues in different
/// symmetry copies.  Note that, ideally, I'd like to test for c3 symmetry, but this is as close as was feasible.
/// @note Can be overriden.
bool
CrosslinkerMoverHelper::selection_is_symmetric(
	core::select::residue_selector::ResidueSubset const & selection,
	core::pose::Pose const &pose
) const {
	//First, get the indices (checking in the process that exactly three residues are selected).
	core::Size r1, r2, r3;
	get_sidechain_indices( selection, r1, r2, r3 );

	//Check whether this is a symmetric pose.
	core::conformation::symmetry::SymmetricConformationCOP conf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
	if ( conf == nullptr ) return false;

	//Get the symmetry information object, which can tell us about relationships between residues:
	core::conformation::symmetry::SymmetryInfoCOP symminfo( conf->Symmetry_Info() );
	if ( symminfo->subunits() != 3 ) return false;

	if ( symminfo->bb_is_independent( r1 ) ) {
		if ( symminfo->bb_follows( r2 ) == r1 && symminfo->bb_follows( r3 ) == r1 ) { return true; }
		else { return false; }
	}
	if ( symminfo->bb_is_independent( r2 ) ) {
		if ( symminfo->bb_follows( r1 ) == r2 && symminfo->bb_follows( r3 ) == r2 ) { return true; }
		else { return false; }
	}
	if ( symminfo->bb_is_independent( r3 ) ) {
		if ( symminfo->bb_follows( r1 ) == r3 && symminfo->bb_follows( r2 ) == r3 ) { return true; }
		else { return false; }
	}
	core::Size const bb_follows_r1( symminfo->bb_follows( r1 ) );
	if ( symminfo->bb_follows( r2 ) == bb_follows_r1 && symminfo->bb_follows( r3 ) == bb_follows_r1 ) { return true; /*Though this shouldn't really be possible.*/ }
	return false;
}

/// @brief Optional steps that the helper can apply before every relaxation round.
/// @details Defaults to doing nothing; can be overriden.  (One example is the TMA helper, which uses this to update amide bond-dependent atom positions).
void
CrosslinkerMoverHelper::pre_relax_round_update_steps(
	core::pose::Pose &/*pose*/,
	core::select::residue_selector::ResidueSubset const &/*selection*/,
	bool const /*whole_structure*/,
	bool const /*symmetric*/,
	bool const /*linker_was_added*/
) const {
	//By default, do nothing.
}

/// @brief Optional steps that the helper can apply after every relaxation round.
/// @details Defaults to doing nothing; can be overriden.  (One example is the TMA helper, which uses this to update amide bond-dependent atom positions).
void
CrosslinkerMoverHelper::post_relax_round_update_steps(
	core::pose::Pose &/*pose*/,
	core::select::residue_selector::ResidueSubset const &/*selection*/,
	bool const /*whole_structure*/,
	bool const /*symmetric*/,
	bool const /*linker_was_added*/
) const {
	//By default, do nothing.
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  Both the symmetric and asymmetric versions call this code.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
CrosslinkerMoverHelper::filter_by_constraints_energy(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	bool const symmetric,
	bool const linker_was_added,
	core::Real const &filter_multiplier
) const {
	core::pose::Pose pose_copy( pose ); //Scoring alters the energies in the pose, so we need a copy here.  We're also going to manipulate it a bit.

	//Remove constraints from pose, and add back only those for the linker:
	protocols::simple_moves::ClearConstraintsMover clear_csts;
	clear_csts.apply(pose_copy);
	if ( symmetric ) {
		add_linker_constraints_symmetric( pose_copy, selection, linker_was_added );
	} else {
		add_linker_constraints_asymmetric( pose_copy, selection );
	}

	core::scoring::ScoreFunctionOP sfxn;
	if ( symmetric ) {
		sfxn = core::scoring::ScoreFunctionOP( core::scoring::symmetry::SymmetricScoreFunctionOP ( new core::scoring::symmetry::SymmetricScoreFunction ) );
	} else {
		sfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
	}

	sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	sfxn->set_weight( core::scoring::angle_constraint, 1.0 );
	sfxn->set_weight( core::scoring::dihedral_constraint, 1.0 );

	(*sfxn)( pose_copy );

	core::Real const cst_energy( pose_copy.energies().total_energy() );
	bool const failure( cst_energy > (50.0 * filter_multiplier) );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Constraint energy for crosslinked pose is " << cst_energy << ".  " << (failure ? "Filter failed." : "Filter passes.") << std::endl;
	}

	return failure;
}



} //crosslinker
} //protocols
} //cyclic_peptide

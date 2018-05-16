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
/// of linkers.
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
#include <protocols/constraint_movers/ClearConstraintsMover.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.CrosslinkerMoverHelper" );

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CrosslinkerMoverHelper::CrosslinkerMoverHelper() :
	symm_type_('A'),
	symm_count_(1)
	//TODO initialize data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
CrosslinkerMoverHelper::CrosslinkerMoverHelper( CrosslinkerMoverHelper const & src ) :
	symm_type_(src.symm_type_),
	symm_count_(src.symm_count_)
	//TODO copy data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
CrosslinkerMoverHelper::~CrosslinkerMoverHelper()= default;

//////////////////////
/// Public Methods ///
//////////////////////

/// @brief Set the symmetry for this crosslinker helper.
void
CrosslinkerMoverHelper::set_symmetry(
	char const symm_type_in,
	core::Size const symm_count_in
) {
	runtime_assert( symm_type_in == 'A' || symm_type_in == 'C' || symm_type_in == 'D' || symm_type_in == 'S' );
	symm_type_ = symm_type_in;
	if ( symm_type_ == 'A' ) {
		runtime_assert( symm_count_in == 1 );
	} else {
		runtime_assert( symm_count_in > 1 );
		if ( symm_type_ == 'S' ) {
			runtime_assert( symm_count_in % 2 == 0 );
		}
	}
	symm_count_ = symm_count_in;
}

/// @brief Given a ResidueSubset with N residues selected, pull out the indices into a vector.
/// @details Overwrites res_indices.
void
CrosslinkerMoverHelper::get_sidechain_indices(
	core::select::residue_selector::ResidueSubset const & selection,
	utility::vector1< core::Size > & res_indices
) const {
	core::Size const nres( selection.size() );
	res_indices.clear();
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( selection[i] ) {
			res_indices.push_back(i);
		}
	}

	runtime_assert_string_msg( res_indices.size(), "Error in protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelper::get_sidechain_indices(): No residues were selected." );
}

/// @brief Determine whether a selection is symmetric.
/// @details Returns true if and only if (a) the pose is symmetric, (b) there are the expected number of symmetry copies, and (c) the selected residues are equivalent residues in different
/// symmetry copies.  Note that, ideally, I'd like to test for CN or SN symmetry, but this is as close as was feasible.
/// @note Can be overriden.
bool
CrosslinkerMoverHelper::selection_is_symmetric(
	core::select::residue_selector::ResidueSubset const & selection,
	core::pose::Pose const &pose,
	core::Size const expected_subunit_count
) const {
	//First, get the indices.
	utility::vector1< core::Size > res_indices;
	get_sidechain_indices( selection, res_indices );

	//Check whether this is a symmetric pose.
	core::conformation::symmetry::SymmetricConformationCOP conf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
	if ( conf == nullptr ) return false;

	//Get the symmetry information object, which can tell us about relationships between residues:
	core::conformation::symmetry::SymmetryInfoCOP symminfo( conf->Symmetry_Info() );
	if ( symminfo->subunits() != expected_subunit_count ) return false;
	if ( (res_indices.size() != expected_subunit_count) && (res_indices.size() % expected_subunit_count != 0) ) return false;

	for ( core::Size i(1), imax(res_indices.size()); i<=imax; ++i ) {
		if ( symminfo->bb_is_independent( res_indices[i] ) ) {
			bool has_symmetry_copy( false );
			for ( core::Size j(1), jmax(res_indices.size()); j<=jmax; ++j ) {
				if ( i==j ) continue;
				if ( symminfo->bb_follows( res_indices[j] ) == res_indices[i] ) {
					has_symmetry_copy = true;
					break;
				}
			}
			if ( !has_symmetry_copy ) return false;
		}
	}

	return true;
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

/// @brief Get the number of expected symmetry subunits, given the symmetry type.
core::Size
CrosslinkerMoverHelper::symm_subunits_expected() const {
	if ( symm_type() == 'D' ) return symm_count()*2;
	return symm_count();
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
	protocols::constraint_movers::ClearConstraintsMover clear_csts;
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

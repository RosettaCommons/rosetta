// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/threefold_linker/ThreefoldLinkerMoverHelper.cc
/// @brief A base class for helper objects that the ThreefoldLinkerMover uses to set up specific types
/// of threefold linkers.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/threefold_linker/ThreefoldLinkerMoverHelper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>

// Protocols headers

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.threefold_linker.ThreefoldLinkerMoverHelper" );

namespace protocols {
namespace cyclic_peptide {
namespace threefold_linker {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ThreefoldLinkerMoverHelper::ThreefoldLinkerMoverHelper() //:
//TODO initialize data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
ThreefoldLinkerMoverHelper::ThreefoldLinkerMoverHelper( ThreefoldLinkerMoverHelper const & /*src*/ ) //:
//TODO copy data here
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
ThreefoldLinkerMoverHelper::~ThreefoldLinkerMoverHelper(){}

//////////////////////
/// Public Methods ///
//////////////////////

/// @brief Given a ResidueSubset with exactly three residues selected, pull out the three indices.
/// @details Overwrites res1, res2, and res3.
void
ThreefoldLinkerMoverHelper::get_sidechain_indices(
	core::select::residue_selector::ResidueSubset const & selection,
	core::Size &res1,
	core::Size &res2,
	core::Size &res3
) const {
	core::Size const nres( selection.size() );
	runtime_assert_string_msg( nres>=3, "Error in protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelper::get_sidechain_indices(): Fewer than three residues are in the pose." );
	res1=0; res2=0; res3=0;
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( selection[i] ) {
			if ( res1==0 ) { res1 = i; }
			else if ( res2==0 ) { res2 = i; }
			else if ( res3==0 ) { res3 = i; }
			else {
				utility_exit_with_message( "Error in protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelper::get_sidechain_indices(): More than three residues were selected." );
			}
		}
	}

	runtime_assert_string_msg( res1>0 && res2>0 && res3>0, "Error in protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelper::get_sidechain_indices(): Fewer than three residues were selected." );
}

/// @brief Determine whether a selection is symmetric.
/// @details Returns true if and only if (a) the pose is symmetric, (b) there are three symmetry copies, and (c) the selected residues are equivalent residues in different
/// symmetry copies.  Note that, ideally, I'd like to test for c3 symmetry, but this is as close as was feasible.
/// @note Can be overriden.
bool
ThreefoldLinkerMoverHelper::selection_is_symmetric(
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


} //threefold_linker
} //protocols
} //cyclic_peptide

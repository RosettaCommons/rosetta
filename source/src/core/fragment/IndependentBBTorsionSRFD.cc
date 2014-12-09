// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/IndependentBBTorsionSRFD.cc
/// @brief  A version of BBTorsionSRFD that considers each torsion independently
///         during is_applicable() and apply() calls when passed a MoveMap (vs
///         the all-torsions-must-be-moveable-or-nothing-is behavior in the
///         original BBTorsionSRFD).
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


// unit headers
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>



namespace core {
namespace fragment {



/// @brief default constructor
IndependentBBTorsionSRFD::IndependentBBTorsionSRFD() :
	Super()
{}


/// @brief constructor
/// @param[in] n_bbtorsions Number of backbone torsions.
/// @param[in] secstruct The single character secondary structure type.
/// @param[in] sequence The single character sequence type.
IndependentBBTorsionSRFD::IndependentBBTorsionSRFD(
	Size const n_bbtorsions,
	char const secstruct,
	char const sequence
) :
	Super( n_bbtorsions, secstruct, sequence )
{}


/// @brief copy constructor
IndependentBBTorsionSRFD::IndependentBBTorsionSRFD( IndependentBBTorsionSRFD const & rval ) :
	Super( rval )
{}


/// @brief default destructor
IndependentBBTorsionSRFD::~IndependentBBTorsionSRFD() {}


/// @brief copy assignment
IndependentBBTorsionSRFD & IndependentBBTorsionSRFD::operator =( IndependentBBTorsionSRFD const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
	}
	return *this;
}


/// @brief clone this object
SingleResidueFragDataOP IndependentBBTorsionSRFD::clone() const {
	return SingleResidueFragDataOP( new IndependentBBTorsionSRFD( *this ) );
}


/// @brief create a new instance of this object
SingleResidueFragDataOP IndependentBBTorsionSRFD::create() const {
	return SingleResidueFragDataOP( new IndependentBBTorsionSRFD() );
}


/// @brief apply only torsions in this fragment marked as moveable in the given
///  MoveMap
/// @param[in] movemap Check for moveable torsions in this MoveMap.
/// @param[in,out] pose The Pose to modify.
/// @param[in] seqpos Insert at this sequence position.
/// @return True if at least one torsion inserted and second level superclass
///  <tt>SecstructSRFD::apply()</tt> succeeded, otherwise false.
/// @details Skips the initial call to primary superclass <tt>BBTorsionSRFD::apply()</tt>
///  because it contains all-or-nothing behavior we want to override.  Instead
///  will go directly to second level <tt>SecstructSRFD::apply()</tt>, and afterwards
///  execute the rest of the procedure.
bool IndependentBBTorsionSRFD::apply(
	MoveMap const & movemap,
	Pose & pose,
	Size const seqpos
) const
{
	// Call SecstructSRFD::apply().  Can't call direct superclass
	// BBTorsionSRFD::apply() because it directly inserts the fragment with
	// no checks, which is something we want to override.
	if ( !Super2::apply( movemap, pose, seqpos ) ) {
		return false; // punt
	}

	// only insert torsions that are allowed to move in the MoveMap
	bool success = false; // at least one torsion inserted?
	for ( Size j = 1, je = nbb(); j <= je; ++j ) {
		id::TorsionID const torsion_id( seqpos, id::BB, j );
		if ( movemap.get( torsion_id ) ) {
			pose.set_torsion( torsion_id, Super::torsion( j ) );
			success = true;
		}
	}

	return success;
}


/// @brief is at least one torsion marked as moveable in the given MoveMap?
/// @param[in] movemap Check for moveable torsions in this MoveMap.
/// @param[in] seqpos Check at this sequence position.
/// @return True if at least one torsion moveable and second level superclass
///  <tt>SecstructSRFD::is_applicable()</tt>, otherwise False.
/// @details Skips the initial call to primary superclass <tt>BBTorsionSRFD::is_applicable()</tt>
///  because it contains all-or-nothing behavior we want to override.  Instead
///  will go directly to second level <tt>SecstructSRFD::is_applicable()</tt>, and afterwards
///  execute the rest of the procedure.
bool IndependentBBTorsionSRFD::is_applicable(
	MoveMap const & movemap,
	Size seqpos
) const
{
	// Call SecstructSRFD::is_applicable().  Can't call direct superclass
	// BBTorsionSRFD::is_applicable() because it checks an all-or-nothing
	// condition we want to override.
	// Only move forward if prior ops successful, otherwise punt.
	if ( !Super2::is_applicable( movemap, seqpos ) ) {
		return false;
	}

	// if any backbone torsion is allowed to move, then this IndependentBBTorsionSRFD
	// is applicable for the given seqpos
	bool applicable = false;
	for ( Size j = 1, je = nbb(); ( !applicable ) && j <= je; ++j ) {
		applicable |= movemap.get( id::TorsionID( seqpos, id::BB, j ) );
	}

	return applicable;
}


} // namespace fragment
} // namespace core

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/GrowLeft.cc
/// @brief instruction to create an n-side extension
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/GrowLeft.hh>

// package headers
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>

// project headers

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief default constructor
GrowLeft::GrowLeft() :
	Super( Interval( 0, 0 ) ),
	pos_( 0 )
{}


/// @brief constructor
/// @param[in] pos grow an n-side extension prior to this position
/// @param[in] ss the secondary structure desired, also defines length of extension
/// @param[in] aa the annotated amino acid sequence desired, default is poly-alanine
/// @param[in] rts the residue type set to use, default FA_STANDARD
/// @remarks length of the *one-letter* aa must equal the length of ss
GrowLeft::GrowLeft(
	Size const pos,
	String const & ss,
	String const & aa,
	ResidueTypeSetCAP rts
) :
	Super( Interval( pos, pos ), rts ),
	pos_( pos ),
	ss_( ss ),
	aa_( aa )
{
	// build poly-alanine if empty string
	if ( aa_.empty() ) {
		aa_ = String( ss_.length(), 'A' );
	}

	debug_assert( ss_.length() == core::pose::annotated_to_oneletter_sequence( aa_ ).length() );
}


/// @brief copy constructor
GrowLeft::GrowLeft( GrowLeft const & rval ) :
	Super( rval ),
	pos_( rval.pos_ ),
	ss_( rval.ss_ ),
	aa_( rval.aa_ )
{}


/// @brief default destructor
GrowLeft::~GrowLeft() {}


/// @brief copy assignment
GrowLeft & GrowLeft::operator =( GrowLeft const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		pos_ = rval.pos_;
		ss_ = rval.ss_;
		aa_ = rval.aa_;
	}
	return *this;
}


/// @brief clone this object
BuildInstructionOP GrowLeft::clone() const {
	return BuildInstructionOP( new GrowLeft( *this ) );
}


/// @brief return a copy of the set of positions within the new region
///  that were pre-existing in the original Pose prior to modify()
/// @return An empty set -- no positions are pre-existing.
GrowLeft::Positions GrowLeft::preexisting_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions that are "new" and did
///  not exist in the original Pose.
/// @return A set of positions spanning the entire region -- all positions
///  are new.
GrowLeft::Positions GrowLeft::new_positions() const {
	using protocols::forge::methods::closed_range;

	Interval const ival = interval();
	return closed_range( ival.left, ival.right );
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has a defined conformation.  E.g. copied residues.
/// @return An empty set -- no positions are defined.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
GrowLeft::Positions GrowLeft::defined_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has an undefined conformation.  E.g. existing or copied residues.
/// @return A set of positions spanning the entire region -- all positions
///  are undefined.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
GrowLeft::Positions GrowLeft::undefined_positions() const {
	// for GrowLeft this is the same as the new positions
	return new_positions();
}


/// @brief return a copy of the MoveMap that defines the moveable/fixed
///  positions/dofs for this instruction
/// @return a MoveMap with [interval.left, interval.right] bb & chi set to true
///  at the MoveMapTorsionID level
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
GrowLeft::MoveMap GrowLeft::movemap() const {
	Positions const newp = new_positions();

	MoveMap mm;

	for ( Positions::const_iterator i = newp.begin(), ie = newp.end(); i != ie; ++i ) {
		mm.set_bb( *i, true );
		mm.set_chi( *i, true );
	}

	return mm;
}


/// @brief update indexing on residue append
void GrowLeft::on_residue_append( LengthEvent const & event ) {
	if ( event.position < pos_ ) {
		//++pos_;
		pos_ += event.length_change;
	}
}


/// @brief update indexing on residue prepend
void GrowLeft::on_residue_prepend( LengthEvent const & event ) {
	if ( event.position <= pos_ ) {
		//++pos_;
		pos_ += event.length_change;
	}
}


/// @brief update indexing on residue delete
void GrowLeft::on_residue_delete( LengthEvent const & event ) {
	if ( pos_ > 1 && event.position <= pos_ && ( int(pos_) + event.length_change > 0) ) {
		//--pos_;
		if ( int(pos_) + event.length_change < int(event.position) ) pos_ = event.position - 1;
		pos_ += event.length_change;
	}
}


/// @brief return the set of positions within the original interval that
///  will be kept in this BuildInstruction
/// @return An empty set -- no positions are kept.
GrowLeft::Positions GrowLeft::original_kept_positions() const {
	return Positions();
}


/// @brief return set of positions within the original interval that will
///  be deleted in this BuildInstruction
/// @return An empty set -- no positions are deleted.
GrowLeft::Positions GrowLeft::original_deleted_positions() const {
	return Positions();
}


/// @brief return set of any fixed positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no fixed positions
GrowLeft::Positions GrowLeft::original_fixed_positions() const {
	Positions fixed;
	fixed.insert( original_interval().right );
	return fixed;
}


/// @brief return set of any mutable positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no mutable positions
GrowLeft::Positions GrowLeft::original_mutable_positions() const {
	return Positions();
}


/// @brief do the actual work of modifying the Pose
void GrowLeft::modify_impl( Pose & pose ) {
	using core::chemical::ResidueTypeCOPs;
	using protocols::forge::methods::grow_left_rtype;
	using protocols::forge::methods::trans_omega;

	// grab residue types from aa string
	ResidueTypeCOPs r_types = core::pose::residue_types_from_sequence( aa_, residue_type_set(), false );
	debug_assert( r_types.size() == ss_.length() );

	// BEGIN POS SHIFT: after this point, pos_ will begin to shift due to length
	// changes in the Pose

	// grow extension
	Size left_endpoint = grow_left_rtype( pose, pos_, r_types.rbegin(), r_types.rend() );
	debug_assert( left_endpoint == pos_ - ss_.length() );

	// END POS SHIFT: after this point, pos_ has stabilized

	// assume proper omega
	Interval new_region = interval();
	trans_omega( new_region.left, new_region.right, pose );

	// set the desired secondary structure
	for ( Size r = left_endpoint, i = 0; r <= pos_ - 1; ++r, ++i ) {
		pose.set_secstruct( r, ss_.at( i ) );
	}
}


/// @brief do the actual reset of intervals, positions, etc to initial state
void GrowLeft::reset_accounting_impl() {
	pos_ = original_interval().right;
}


} // namespace build
} // namespace forge
} // namespace protocols

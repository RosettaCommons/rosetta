// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/build/GrowRight.cc
/// @brief instruction to create a c-side extension
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/GrowRight.hh>

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
GrowRight::GrowRight() :
	Super( Interval( 0, 0 ) ),
	pos_( 0 )
{}


/// @brief constructor
/// @param[in] pos grow a c-side extension after to this position
/// @param[in] ss the secondary structure desired, also defines length of extension
/// @param[in] aa the annotated amino acid sequence desired, default is poly-alanine
/// @param[in] rts the residue type set to use, default FA_STANDARD
/// @remarks length of the *one-letter* aa must equal the length of ss
GrowRight::GrowRight(
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

	assert( ss_.length() == core::pose::annotated_to_oneletter_sequence( aa_ ).length() );
}


/// @brief copy constructor
GrowRight::GrowRight( GrowRight const & rval ) :
	Super( rval ),
	pos_( rval.pos_ ),
	ss_( rval.ss_ ),
	aa_( rval.aa_ )
{}


/// @brief default destructor
GrowRight::~GrowRight() {}


/// @brief copy assignment
GrowRight & GrowRight::operator =( GrowRight const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		pos_ = rval.pos_;
		ss_ = rval.ss_;
		aa_ = rval.aa_;
	}
	return *this;
}


/// @brief clone this object
BuildInstructionOP GrowRight::clone() const {
	return BuildInstructionOP( new GrowRight( *this ) );
}


/// @brief return a copy of the set of positions within the new region
///  that were pre-existing in the original Pose prior to modify()
/// @return An empty set -- no positions are pre-existing.
GrowRight::Positions GrowRight::preexisting_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions that are "new" and did
///  not exist in the original Pose.
/// @return A set of positions spanning the entire region -- all positions
///  are new.
GrowRight::Positions GrowRight::new_positions() const {
	using protocols::forge::methods::closed_range;

	Interval const ival = interval();
	return closed_range( ival.left, ival.right );
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has a defined conformation.  E.g. existing or copied residues.
/// @return An empty set -- no positions are defined.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
GrowRight::Positions GrowRight::defined_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has an undefined conformation.  E.g. newly created residues.
/// @return A set of positions spanning the entire region -- all positions
///  are undefined.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
GrowRight::Positions GrowRight::undefined_positions() const {
	// for GrowRight this is the same as the new positions
	return new_positions();
}


/// @brief return a copy of the MoveMap that defines the moveable/fixed
///  positions/dofs for this instruction
/// @return a MoveMap with [interval.left, interval.right] bb & chi set to true
///  at the MoveMapTorsionID level
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
GrowRight::MoveMap GrowRight::movemap() const {
	Positions const newp = new_positions();

	MoveMap mm;

	for ( Positions::const_iterator i = newp.begin(), ie = newp.end(); i != ie; ++i ) {
		mm.set_bb( *i, true );
		mm.set_chi( *i, true );
	}

	return mm;
}


/// @brief update indexing on residue append
void GrowRight::on_residue_append( LengthEvent const & event ) {
	if ( event.position < pos_ ) {
		//++pos_;
		pos_ += event.length_change;
	}
}


/// @brief update indexing on residue prepend
void GrowRight::on_residue_prepend( LengthEvent const & event ) {
	if ( event.position <= pos_ ) {
		//++pos_;
		pos_ += event.length_change;
	}
}


/// @brief update indexing on residue delete
void GrowRight::on_residue_delete( LengthEvent const & event ) {
	if ( event.position <= pos_ ) {
		//--pos_;
		if ( int(pos_) + event.length_change < int(event.position) ) pos_ = event.position;
		pos_ += event.length_change;
	}
}


/// @brief return the set of positions within the original interval that
///  will be kept in this BuildInstruction
/// @return An empty set -- no positions are kept.
GrowRight::Positions GrowRight::original_kept_positions() const {
	return Positions();
}


/// @brief return set of positions within the original interval that will
///  be deleted in this BuildInstruction
/// @return An empty set -- no positions are deleted.
GrowRight::Positions GrowRight::original_deleted_positions() const {
	return Positions();
}


/// @brief return set of any fixed positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no fixed positions necessary
GrowRight::Positions GrowRight::original_fixed_positions() const {
	Positions fixed;
	fixed.insert( original_interval().left );
	return fixed;
}


/// @brief return set of any mutable positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no mutable positions
GrowRight::Positions GrowRight::original_mutable_positions() const {
	return Positions();
}


/// @brief do the actual work of modifying the Pose
void GrowRight::modify_impl( Pose & pose ) {
	using core::chemical::ResidueTypeCOPs;
	using protocols::forge::methods::grow_right_rtype;
	using protocols::forge::methods::trans_omega;

	// grab residue types from aa string
	ResidueTypeCOPs r_types = core::pose::residue_types_from_sequence( aa_, residue_type_set(), false );
	assert( r_types.size() == ss_.length() );

	// BEGIN POS SHIFT: after this point, pos_ will begin to shift due to length
	// changes in the Pose

	// grow extension
	Size right_endpoint = grow_right_rtype( pose, pos_, r_types.begin(), r_types.end() );
	assert( right_endpoint == pos_ + ss_.length() );

	// END POS SHIFT: after this point, pos_ has stabilized

	// assume proper omega
	Interval new_region = interval();
	Size const omega_left = ( new_region.left == 1 ) ? new_region.left : new_region.left - 1;
	trans_omega( omega_left, new_region.right, pose );

	// set the desired secondary structure
	for ( Size r = pos_ + 1, i = 0; r <= right_endpoint; ++r, ++i ) {
		pose.set_secstruct( r, ss_.at( i ) );
	}
}


/// @brief do the actual reset of intervals, positions, etc to initial state
void GrowRight::reset_accounting_impl() {
	pos_ = original_interval().left;
}


} // namespace build
} // namespace forge
} // namespace protocols

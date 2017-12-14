// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/build/Bridge.cc
/// @brief  connect two contiguous but disjoint sections of a Pose into one
///         continuous section
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/Bridge.hh>

// package headers
#include <protocols/forge/build/GrowRight.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>

// project headers

#include <core/conformation/util.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

#include <core/pose/annotated_sequence.hh>
#include <utility>
#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief default constructor
Bridge::Bridge() :
	Super()
{}


/// @brief sec.struct only constructor (poly-alanine for new region)
/// @param[in] interval bridge these two residues
/// @param[in] ss the secondary structure desired, also defines length of new bridge,
///  region between the anchor positions, can be empty
/// @param[in] rts the residue type set to use, default FA_STANDARD
/// @remarks length of the *one-letter* aa must equal the length of ss
Bridge::Bridge(
	Interval const & i,
	String const & ss,
	ResidueTypeSetCAP rts
) :
	Super( i, rts ),
	interval_( i ),
	ss_( ss )
{
	// build poly-alanine if empty string
	if ( !ss_.empty() && aa_.empty() ) {
		aa_ = String( ss_.length(), 'A' );
	}

	debug_assert( ss_.length() == core::pose::annotated_to_oneletter_sequence( aa_ ).length() );
}


/// @brief full constructor
/// @param[in] interval bridge these two residues
/// @param[in] ss the secondary structure desired, also defines length of new bridge,
///  region between the anchor positions, can be empty
/// @param[in] aa the annotated amino acid sequence desired, default is poly-alanine
/// @param[in] rts the residue type set to use, default FA_STANDARD
/// @remarks length of the *one-letter* aa must equal the length of ss
Bridge::Bridge(
	Interval const & i,
	String const & ss,
	String const & aa,
	ResidueTypeSetCAP rts
) :
	Super( i, rts ),
	interval_( i ),
	ss_( ss ),
	aa_( aa )
{
	// build poly-alanine if empty string
	if ( !ss_.empty() && aa_.empty() ) {
		aa_ = String( ss_.length(), 'A' );
	}

	debug_assert( ss_.length() == core::pose::annotated_to_oneletter_sequence( aa_ ).length() );
}


/// @brief copy constructor
Bridge::Bridge( Bridge const & /*rval*/ ) = default;


/// @brief default destructor
Bridge::~Bridge() = default;


/// @brief copy assignment
Bridge & Bridge::operator =( Bridge const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		interval_ = rval.interval_;
		ss_ = rval.ss_;
		aa_ = rval.aa_;
	}
	return *this;
}


/// @brief clone this object
BuildInstructionOP Bridge::clone() const {
	return BuildInstructionOP( new Bridge( *this ) );
}


/// @brief return a copy of the set of positions within the new region
///  that were pre-existing in the original Pose prior to modify()
/// @return A set containing two positions -- interval.left and interval.right.
Bridge::Positions Bridge::preexisting_positions() const {
	Positions pre;

	// the interval endpoints are defined
	Interval const ival = interval();
	pre.insert( ival.left );
	pre.insert( ival.right );

	return pre;
}


/// @brief return a copy of the set of positions that are "new" and did
///  not exist in the original Pose.
/// @return A set containing positions spanning [interval.left+1, interval.right-1].
Bridge::Positions Bridge::new_positions() const {
	using protocols::forge::methods::closed_range;

	// everything except for the interval endpoints is undefined
	Interval const ival = interval();
	return closed_range( ival.left + 1, ival.right - 1 );
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has a defined conformation.  E.g. existing or copied residues.
/// @return A set containing two positions -- interval.left and interval.right.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
Bridge::Positions Bridge::defined_positions() const {
	// for Bridge this is the same as the pre-existing positions
	return preexisting_positions();
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has an undefined conformation.  E.g. newly created residues.
/// @return A set containing positions spanning [interval.left+1, interval.right-1].
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
Bridge::Positions Bridge::undefined_positions() const {
	// for Bridge this is the same as the new positions
	return new_positions();
}


/// @brief return a copy of the MoveMap that defines the moveable/fixed
///  positions/dofs for this instruction
/// @return a MoveMap with [interval.left, interval.right] bb & chi set to true
///  at the MoveMapTorsionID level
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
Bridge::MoveMap Bridge::movemap() const {
	// the entire region is mutable
	Interval const ival = interval();

	MoveMap mm;

	for ( Size i = ival.left; i <= ival.right; ++i ) {
		mm.set_bb( i, true );
		mm.set_chi( i, true );
	}

	return mm;
}


/// @brief update indexing on residue append
/// @remarks left and right endpoints of the interval can travel independently
void Bridge::on_residue_append( LengthEvent const & event ) {
	if ( event.position < interval_.left ) {
		//++interval_.left;
		interval_.left = interval_.left + event.length_change;
	}

	if ( event.position < interval_.right ) {
		//++interval_.right;
		interval_.right +=  event.length_change;
	}
}


/// @brief update indexing on residue prepend
/// @remarks left and right endpoints of the interval can travel independently
void Bridge::on_residue_prepend( LengthEvent const & event ) {
	if ( event.position <= interval_.left ) {
		//++interval_.left;
		interval_.left += event.length_change;
	}

	if ( event.position <= interval_.right ) {
		//++interval_.right;
		interval_.right +=  event.length_change;
	}
}


/// @brief update indexing on residue delete
/// @remarks Left and right endpoints of the interval can travel independently.
void Bridge::on_residue_delete( LengthEvent const & event ) {
	// event.position == interval.left is not caught below.
	// It has context dependent consequences and is manually corrected for
	// during modify().
	if ( event.position < interval_.left ) { // left
		//--interval_.left;
		if ( int(interval_.left) + event.length_change < int(event.position) ) interval_.left = event.position;
		else interval_.left += event.length_change;
	}

	if ( event.position < interval_.right ) { // right
		//--interval_.right;
		if ( int(interval_.right) + event.length_change < int(event.position) ) interval_.right = event.position;
		else interval_.right += event.length_change;
	}
}


/// @brief return the set of positions within the original interval that
///  will be kept in this BuildInstruction
/// @return A set containing the endpoints of the original interval.
Bridge::Positions Bridge::original_kept_positions() const {
	Positions kept;
	kept.insert( original_interval().left );
	kept.insert( original_interval().right );
	return kept;
}


/// @brief return set of positions within the original interval that will
///  be deleted in this BuildInstruction
/// @return A set containing the positions in [original.left+1, original.right-1].
Bridge::Positions Bridge::original_deleted_positions() const {
	using protocols::forge::methods::closed_range;
	return closed_range( original_interval().left + 1, original_interval().right - 1 );
}


/// @brief return set of any fixed positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no fixed positions necessary
Bridge::Positions Bridge::original_fixed_positions() const {
	Positions fixed;
	fixed.insert( original_interval().left - 1 );
	fixed.insert( original_interval().right + 1 );
	return fixed;
}


/// @brief return set of any mutable positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no mutable positions
Bridge::Positions Bridge::original_mutable_positions() const {
	Positions muta;
	muta.insert( original_interval().left );
	muta.insert( original_interval().right );
	return muta;
}


/// @brief do the actual work of modifying the Pose
void Bridge::modify_impl( Pose & pose ) {
	using core::kinematics::FoldTree;
	using protocols::forge::build::GrowRight;

	using core::pose::remove_lower_terminus_type_from_pose_residue;
	using core::pose::remove_upper_terminus_type_from_pose_residue;
	using core::conformation::idealize_position;
	using protocols::forge::methods::remove_cutpoint;
	using protocols::forge::methods::remove_cutpoint_variants;

	// check conditions
	runtime_assert( interval_.left > 0 );
	runtime_assert( interval_.right <= pose.size() );
	runtime_assert( interval_.left == interval_.right - 1 );
	runtime_assert( pose.fold_tree().is_cutpoint( interval_.left ) );

	// safety, remove any cutpoint variants in case they exist
	remove_cutpoint_variants( pose, interval_.left );
	remove_cutpoint_variants( pose, interval_.right );

	// safety, remove any termini variants in case they exist
	if ( pose.residue( interval_.left ).is_upper_terminus() ) {
		core::pose::remove_upper_terminus_type_from_pose_residue( pose, interval_.left );
	}

	if ( pose.residue( interval_.right ).is_lower_terminus() ) {
		core::pose::remove_lower_terminus_type_from_pose_residue( pose, interval_.right );
	}

	// BEGIN INTERVAL SHIFT: after this point, interval_ will begin to shift due to length
	// changes in the Pose

	if ( !ss_.empty() ) {
		// Grow out the new section from interval_.left using GrowRight.
		// GrowRight will also handle setting proper omega and secondary
		// structure.
		GrowRight grow( interval_.left, ss_, aa_, residue_type_set().get_self_ptr() );
		grow.modify( pose );
	}

	// END INTERVAL SHIFT: after this point, interval_ has stabilized and stores
	// the new endpoints of the rebuilt segment

	// seal the tree at the cutpoint between [left, right]
	FoldTree new_ft = pose.fold_tree();
	remove_cutpoint( interval_.right - 1, new_ft );
	pose.fold_tree( new_ft );

	// the anchor positions will need to move during remodeling (e.g. fragment
	// insertion), so go ahead and idealize them
	idealize_position( interval_.left, pose.conformation() );
	idealize_position( interval_.right, pose.conformation() );

	// assume proper omega
	pose.set_omega( interval_.left, 180.0 );
	pose.set_omega( interval_.right, 180.0 );
}


/// @brief do the actual reset of intervals, positions, etc to initial state
void Bridge::reset_accounting_impl() {
	interval_ = original_interval();
}


} // namespace build
} // namespace forge
} // namespace protocols

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/build/SegmentSwap.cc
/// @brief instruction to swap a segment with an external segment
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/SegmentSwap.hh>

// package headers
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/methods/util.hh>

// project headers

#include <core/conformation/Conformation.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>

#include <core/conformation/util.hh>
#include <utility/vector1.hh>



namespace protocols {
namespace forge {
namespace build {


/// @brief default constructor
SegmentSwap::SegmentSwap() :
	Super( Interval( 0.0, 0.0 ) ),
	interval_( 0.0, 0.0 )
{}


/// @brief constructor
/// @param[in] interval swap out this range of residues
/// @param[in] move_map fixed backbone residues in this movemap will be used for new jumps
/// @param[in] swap_in swap in this pose
/// @remarks Procedure will attempt to honor the movemap as much as it can.
///  The caveat is that sequences of calls to some FoldTree routines may shift
///  the jumps internally in a way that is not easily predictable.  If the
///  procedure cannot  find an allowed residue for a jump, it will make a jump
///  to the (lower) median residue in the disconnected fold tree interval.
SegmentSwap::SegmentSwap(
	Interval const & i,
	MoveMap const & swap_in_movemap,
	Pose const & swap_in
) :
	Super( i, swap_in.residue( 1 ).residue_type_set().get_self_ptr() ),
	interval_( i ),
	swap_in_movemap_( swap_in_movemap ),
	swap_in_( swap_in )
{
	init();
}


/// @brief copy constructor
SegmentSwap::SegmentSwap( SegmentSwap const & rval ) :
	Super( rval ),
	interval_( rval.interval_ ),
	swap_in_movemap_( rval.swap_in_movemap_ ),
	swap_in_( rval.swap_in_ )
{}


/// @brief default destructor
SegmentSwap::~SegmentSwap() {}


/// @brief copy assignment
SegmentSwap & SegmentSwap::operator =( SegmentSwap const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		interval_ = rval.interval_;
		swap_in_movemap_ = rval.swap_in_movemap_;
		swap_in_ = rval.swap_in_;
	}
	return *this;
}


/// @brief clone this object
BuildInstructionOP SegmentSwap::clone() const {
	return new SegmentSwap( *this );
}


/// @brief fixed backbone residues in this movemap will be used for new jumps
/// @remarks Procedure will attempt to honor this movemap as much as it can.
///  The caveat is that sequences of calls to some FoldTree routines may shift
///  the jumps internally in a way that is not easily predictable.  If the
///  procedure cannot  find an allowed residue for a jump, it will make a jump
///  to the (lower) median residue in the disconnected fold tree interval.
SegmentSwap::MoveMap const & SegmentSwap::swap_in_movemap() const {
	return swap_in_movemap_;
}


/// @brief the pose to swap in
SegmentSwap::Pose const & SegmentSwap::swap_in() const {
	return swap_in_;
}


/// @brief a copy of the working range of residues specifying the swapped region
/// @details This residue range can change wrt length changes in Pose /Conformation
///  being watched.
Interval SegmentSwap::interval() const {
	return interval_;
}


/// @brief return a copy of the set of positions within the new region
///  that were pre-existing in the original Pose prior to modify()
/// @return An empty set -- no positions are pre-existing.
SegmentSwap::Positions SegmentSwap::preexisting_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions that are "new" and did
///  not exist in the original Pose.
/// @return A set of positions spanning the interval -- all positions are
///  are defined.
SegmentSwap::Positions SegmentSwap::new_positions() const {
	using protocols::forge::methods::closed_range;

	Interval const ival = interval();
	return closed_range( ival.left, ival.right );
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has a defined conformation.  E.g. existing or copied residues.
/// @return A set of positions spanning the interval -- all positions are
///  are defined.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentSwap::Positions SegmentSwap::defined_positions() const {
	// for SegmentSwap this is the same as the new positions
	return new_positions();
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has an undefined conformation.  E.g. newly created residues.
/// @return An empty set -- no undefined positions.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentSwap::Positions SegmentSwap::undefined_positions() const {
	return Positions();
}


/// @brief return a copy of the MoveMap that defines the moveable/fixed
///  positions/dofs for this instruction
/// @return a MoveMap with [interval.left, interval.right] bb set to false
///  at the MoveMapTorsionID level
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentSwap::MoveMap SegmentSwap::movemap() const {
	MoveMap mm;

	Positions const newp = new_positions();
	for ( Positions::const_iterator i = newp.begin(), ie = newp.end(); i != ie; ++i ) {
		mm.set_bb( *i, false );
	}

	return mm;
}


/// @brief update indexing on residue append
void SegmentSwap::on_residue_append( LengthEvent const & event ) {
	if ( event.position < interval_.left ) {
		//++interval_.left;
		interval_.left += event.length_change;
	}

	if ( event.position < interval_.right ) {
		++interval_.right;
		//interval_.right += event.length_change;
	}
}


/// @brief update indexing on residue prepend
void SegmentSwap::on_residue_prepend( LengthEvent const & event ) {
	if ( event.position <= interval_.left ) {
		//++interval_.left;
		interval_.left += event.length_change;
	}

	if ( event.position <= interval_.right ) {
		//++interval_.right;
		interval_.right += event.length_change;
	}
}


/// @brief update indexing on residue delete
void SegmentSwap::on_residue_delete( LengthEvent const & event ) {
	if ( event.position < interval_.left ) { // left
		//--interval_.left;
		if( int(interval_.left) + event.length_change < int(event.position) ) interval_.left = event.position;
		else interval_.left += event.length_change;
	}

	if ( event.position < interval_.right ) { // right
		//--interval_.right;
		if( int(interval_.right) + event.length_change < int(event.position) ) interval_.right = event.position;
		else interval_.right += event.length_change;
	}
}


/// @brief return the set of positions within the original interval that
///  will be kept in this BuildInstruction
/// @return An empty set -- no positions are kept.
SegmentSwap::Positions SegmentSwap::original_kept_positions() const {
	return Positions();
}


/// @brief return set of positions within the original interval that will
///  be deleted in this BuildInstruction
/// @return An empty set -- no positions are deleted.
SegmentSwap::Positions SegmentSwap::original_deleted_positions() const {
	using protocols::forge::methods::closed_range;
	return closed_range( original_interval().left, original_interval().right );
	//	return Positions();
}


/// @brief return set of any fixed positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no fixed positions necessary
SegmentSwap::Positions SegmentSwap::original_fixed_positions() const {
	Positions fixed;
	return fixed;
}


/// @brief return set of any mutable positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
SegmentSwap::Positions SegmentSwap::original_mutable_positions() const {
	using protocols::forge::methods::closed_range;
	return closed_range( original_interval().left, original_interval().right );
}


/// @brief do the actual work of modifying the Pose
void SegmentSwap::modify_impl( Pose & pose ) {
	using core::kinematics::FoldTree;
	using core::pose::PDBInfo;
	using core::pose::PDBInfoOP;
	using protocols::forge::methods::replace;

	// save interim FoldTree for passing into replace() later
	FoldTree old_ft = pose.fold_tree();

	// BEGIN INTERVAL SHIFT: after this point, interval_ will begin to shift due to length
	// changes in the Pose

	// Blow away the old segment.  delete_residue_range_slow() should keep
	// all coordinates intact, though the bond lengths/angles adjacent to
	// the deletion might get funky.
	pose.conformation().delete_residue_range_slow( interval_.left, interval_.right );
	assert( interval_.left == interval_.right );

	assert( interval_.left > 1 ); // conformation does not appear to handle certain modifications at 1

	// set temporary fold tree with cut introduced so we can insert residues
	// by jump
	FoldTree tmp_ft = pose.fold_tree();
	tmp_ft.new_jump( interval_.left - 1, interval_.left, interval_.left - 1 );
	pose.fold_tree( tmp_ft );

	// temporarily shift interval_.left back by one so it sits just prior
	// to modification region and append modifications below won't affect it
	--interval_.left;

	// Swap the entire swap_in_ Pose by jump.
	Size append_pos = interval_.left + 1;
	for ( Size i = 1, ie = swap_in_.n_residue(); i <= ie; ++i ) {
		// temporarily add via jump from interval_.left
		pose.insert_residue_by_jump( swap_in_.residue( i ), append_pos, interval_.left, "", "" );
		++append_pos;
	}

	assert( append_pos == interval_.right );
	assert( append_pos == interval_.left + swap_in_.n_residue() + 1 );

	// correct interval_;
	++interval_.left;
	--interval_.right;

	// END INTERVAL SHIFT: after this point, interval_ has stabilized and stores
	// the new endpoints of the rebuilt segment

	assert( ( interval_.right - interval_.left + 1 ) == swap_in_.n_residue() );

	// construct fold tree
	FoldTree new_ft = replace( old_ft, original_interval().left, original_interval().right,
	                           swap_in_movemap_, swap_in_.fold_tree() );

	// set proper topology
	pose.fold_tree( new_ft );

	// copy secstruct
	for ( Size i = interval_.left, j = 1; i <= interval_.right; ++i, ++j ) {
		pose.set_secstruct( i, swap_in_.secstruct( j ) );
	}

	// copy pdb info if it exists; PDBInfo should always be
	// obsolete after leaving modify_impl()
	if ( swap_in_.pdb_info().get() != NULL ) {

		if ( pose.pdb_info().get() == NULL ) {
			pose.pdb_info( new PDBInfo( pose ) );
		}

		// force obsolete
		pose.pdb_info()->obsolete( true );

		// copy
		pose.pdb_info()->copy( *swap_in_.pdb_info(), 1, swap_in_.n_residue(), interval_.left );
	}
}


/// @brief do the actual reset of intervals, positions, etc to initial state
void SegmentSwap::reset_accounting_impl() {
	interval_ = original_interval();
}


/// @brief init to be called during non-default constructors
void SegmentSwap::init() {
	using core::conformation::remove_lower_terminus_type_from_conformation_residue;
	using core::conformation::remove_upper_terminus_type_from_conformation_residue;

	// remove lower/upper terminus only at 1, nres
	if ( swap_in_.n_residue() > 0 ) {
		if ( swap_in_.residue( 1 ).is_lower_terminus() ) {
			core::conformation::remove_lower_terminus_type_from_conformation_residue( swap_in_.conformation(), 1 );
		}

		if ( swap_in_.residue( swap_in_.n_residue() ).is_upper_terminus() ) {
			core::conformation::remove_upper_terminus_type_from_conformation_residue( swap_in_.conformation(), swap_in_.n_residue() );
		}
	}
}



} // namespace build
} // namespace forge
} // namespace protocols

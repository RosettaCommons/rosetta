// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/ConnectRight.cc
/// @brief instruction to connect one Pose onto the right side of another
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/ConnectRight.hh>

// package headers
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/methods/util.hh>

// rosetta headers
#include <core/id/AtomID.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

// C++ headers
#include <limits>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief default constructor
ConnectRight::ConnectRight() :
	Super( Interval( 0, 0 ) ),
	left_position_( 0 ),
	right_position_( 0 ),
	use_rt_( false ),
	interval_( 0.0, 0.0 )
{
	init();
}


/// @brief position to position jump constructor
/// @param[in] left_position connect at this position on 'pose_left'
///  passed into modify()
/// @param[in] right_position connect at this position on 'pose_right'
/// @param[in] pose_right connect this pose to the right of pose_left when
///  modify( pose_left ) is called
/// @param[in] start_new_chain start new chain when connecting?  does not
///  handle termini variants
ConnectRight::ConnectRight(
	Size const left_position,
	Size const right_position,
	Pose const & pose_right
) :
	Super( Interval( left_position, left_position ), pose_right.residue_type_set_for_pose() ),
	left_position_( left_position ),
	right_position_( right_position ),
	pose_right_( pose_right ),
	use_rt_( false ),
	interval_( 1, pose_right.size() )
{
	init();
}


/// @brief copy constructor
ConnectRight::ConnectRight( ConnectRight const & rval ) :
	Super( rval ),
	left_position_( rval.left_position_ ),
	right_position_( rval.right_position_ ),
	pose_right_( rval.pose_right_ ),
	use_rt_( rval.use_rt_ ),
	left_stub_atoms_( rval.left_stub_atoms_ ),
	right_stub_atoms_( rval.right_stub_atoms_ ),
	rt_( rval.rt_ ),
	interval_( rval.interval_ )
{}


/// @brief default destructor
ConnectRight::~ConnectRight() {}


/// @brief copy assignment
ConnectRight & ConnectRight::operator =( ConnectRight const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		left_position_ = rval.left_position_;
		right_position_ = rval.right_position_;
		pose_right_ = rval.pose_right_;
		use_rt_ = rval.use_rt_;
		left_stub_atoms_ = rval.left_stub_atoms_;
		right_stub_atoms_ = rval.right_stub_atoms_;
		rt_ = rval.rt_;
		interval_ = rval.interval_;
	}
	return *this;
}


/// @brief clone this object
BuildInstructionOP ConnectRight::clone() const {
	return BuildInstructionOP( new ConnectRight( *this ) );
}


/// @brief connect this pose to the right of pose_left when modify( pose_left )
///  is called
ConnectRight::Pose const & ConnectRight::pose_right() {
	return pose_right_;
}


/// @brief return a copy of the set of positions within the new region
///  that were pre-existing in the original Pose prior to modify()
/// @return An empty set -- no positions are pre-existing.
ConnectRight::Positions ConnectRight::preexisting_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions that are "new" and did
///  not exist in the original Pose.
/// @return A set spanning the entire interval -- all positions are new.
ConnectRight::Positions ConnectRight::new_positions() const {
	using protocols::forge::methods::closed_range;

	Interval const ival = interval();
	return closed_range( ival.left, ival.right );
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has a defined conformation.  E.g. existing or copied residues.
/// @return A set spanning the entire interval -- all positions are defined.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
ConnectRight::Positions ConnectRight::defined_positions() const {
	// for ConnectRight this is the same as the new positions
	return new_positions();
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has an undefined conformation.  E.g. newly created residues.
/// @return An empty set -- no undefined positions.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
ConnectRight::Positions ConnectRight::undefined_positions() const {
	return Positions();
}


/// @brief return a copy of the MoveMap that defines the moveable/fixed
///  positions/dofs for this instruction
/// @return a MoveMap with [interval.left, interval.right] bb set to false
///  at the MoveMapTorsionID level
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
ConnectRight::MoveMap ConnectRight::movemap() const {
	MoveMap mm;

	Positions const newp = new_positions();
	for ( Positions::const_iterator i = newp.begin(), ie = newp.end(); i != ie; ++i ) {
		mm.set_bb( *i, false );
	}

	return mm;
}


/// @brief extract appropriately computed transform between two stubs that
///  represent a jump in the given Pose between the two residues and set it
///  as the rt for this ConnectRight
/// @param[in] pose The Pose to use.
/// @param[in] jump_start_residue The starting residue of the jump.
///  This residue should be the equivalent of the jump position on pose_left.
/// @param[in] jump_stop-residue The stopping residue of the jump.
///  This residue should be the equivalent of the jump position on pose_right.
/// @remarks Uses left_stub_atoms() and right_stub_atoms() as the stub atoms
///  to compute the transform.  Remember to set use_rt() to True after
///  calling this function if you want to actually use the transform
///  during modify().
void ConnectRight::extract_rt(
	Pose const & pose,
	Size const jump_start_residue,
	Size const jump_stop_residue
)
{
	using core::id::StubID;

	assert( jump_start_residue != jump_stop_residue );
	assert( jump_start_residue <= pose.size() );
	assert( jump_stop_residue <= pose.size() );

	StubID  left_stub_id( core::pose::named_stub_id_to_stub_id(  left_named_stub_id( jump_start_residue ), pose ) );
	StubID right_stub_id( core::pose::named_stub_id_to_stub_id( right_named_stub_id( jump_stop_residue  ), pose ) );

	rt_ = pose.conformation().get_stub_transform( left_stub_id, right_stub_id );
}


/// @brief update indexing on residue append
void ConnectRight::on_residue_append( LengthEvent const & event ) {
	if ( event.position < left_position_ ) {
		//++left_position_;
		left_position_ = left_position_ + event.length_change;
	}

	if ( modify_was_successful() && event.position < interval_.left ) {
		//++interval_.left;
		//++interval_.right;
		interval_.left = interval_.left + event.length_change;
		interval_.right = interval_.right + event.length_change;
	}
}


/// @brief update indexing on residue prepend
void ConnectRight::on_residue_prepend( LengthEvent const & event ) {
	if ( event.position <= left_position_ ) {
		//++left_position_;
		left_position_ = left_position_ + event.length_change;
	}

	if ( modify_was_successful() && event.position <= interval_.left ) {
		//++interval_.left;
		//++interval_.right;
		interval_.left = interval_.left + event.length_change;
		interval_.right = interval_.right + event.length_change;

	}
}


/// @brief update indexing on residue delete
void ConnectRight::on_residue_delete( LengthEvent const & event ) {
	assert( event.position != left_position_ );
	if ( event.position < left_position_ ) {
		//--left_position_;
		if ( int(left_position_) + event.length_change < int(event.position) ) left_position_ = event.position;
		else left_position_ = left_position_ + event.length_change;
	}

	if ( modify_was_successful() && event.position < interval_.left ) {
		//--interval_.left;
		//--interval_.right;
		if ( int(interval_.left) + event.length_change < int(event.position) ) {
			interval_.right = interval_.right - (interval_.left - event.position);
			interval_.left = event.position;
		} else {
			interval_.left = interval_.left + event.length_change;
			interval_.right = interval_.right + event.length_change;
		}
	}
}


/// @brief return the set of positions within the original interval that
///  will be kept in this BuildInstruction
/// @return An empty set -- no positions are kept.
ConnectRight::Positions ConnectRight::original_kept_positions() const {
	return Positions();
}


/// @brief return set of positions within the original interval that will
///  be deleted in this BuildInstruction
/// @return An empty set -- no positions are deleted.
ConnectRight::Positions ConnectRight::original_deleted_positions() const {
	return Positions();
}


/// @brief return set of any fixed positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no fixed positions
ConnectRight::Positions ConnectRight::original_fixed_positions() const {
	Positions fixed;
	fixed.insert( original_interval().left );
	return fixed;
}


/// @brief return set of any mutable positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no mutable positions
ConnectRight::Positions ConnectRight::original_mutable_positions() const {
	return Positions();
}


/// @brief do the actual work of modifying the Pose
void ConnectRight::modify_impl( Pose & pose_left ) {
	using core::id::StubID;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;
	using core::pose::PDBInfo;
	using core::pose::PDBInfoOP;
	using protocols::forge::methods::merge;

	// cache data
	Size const original_left_nres = pose_left.size();

	// construct fold tree
	FoldTree new_ft = merge( pose_left.fold_tree(), left_position_, pose_right_.fold_tree(), right_position_ );

	// BEGIN POSE MODIFY: after this point conformation of pose_left is modified

	// Create [pose_left, pose_right] by attaching all pose_right residues via jump.
	core::Size current_chain = std::numeric_limits< int >::min();

	// if the last residue of pose_left is not an upper terminus and the first
	// residue of pose_right is not a lower terminus, then assume that we're
	// not starting a new chain
	if (
			!pose_left.residue( original_left_nres ).is_upper_terminus() &&
			!pose_right_.residue( 1 ).is_lower_terminus()
			) {
		current_chain = pose_right_.residue( 1 ).chain();
	}

	// Doing the copy with append_residue_by_jump currently appears to be pretty
	// expensive.  However, it's probably the safest way since we might be
	// dealing with Poses that might be "non-standard/weird" vs regular
	// polymeric (protein/dna/rna) Poses.
	for ( Size i = 1, ie = pose_right_.size(); i <= ie; ++i ) {
		if ( pose_right_.residue( i ).chain() != current_chain ) {
			current_chain = pose_right_.residue( i ).chain();
			pose_left.append_residue_by_jump( pose_right_.residue( i ), pose_left.size(), "", "", true );
		} else if ( i > 1 && pose_right_.residue( i ).is_polymer_bonded( i-1 ) ) {
			//If it's bonded in the source pose we probably want it bonded in the destination pose
			pose_left.append_residue_by_bond( pose_right_.residue( i ) );
		} else {
			pose_left.append_residue_by_jump( pose_right_.residue( i ), pose_left.size(), "", "", false );
		}
	}

	// set correct topology
	pose_left.fold_tree( new_ft );

	// MIDDLE POSE MODIFY: after this point non-conformation info of pose_left
	// is modified

	// transfer sec.struct
	for ( Size i = 1, ie = pose_right_.size(); i <= ie; ++i ) {
		pose_left.set_secstruct( i + original_left_nres, pose_right_.secstruct( i ) );
	}

	// Handle PDBInfo separately. PDBInfo in pose_left
	// should always be obsolete after leaving modify_impl()
	// as a safety.
	if ( pose_right_.pdb_info().get() != NULL ) {

		// if pose_left doesn't have PDBInfo, create it
		if ( pose_left.pdb_info().get() == NULL ) {
			pose_left.pdb_info( PDBInfoOP( new PDBInfo( pose_left ) ) );
		}

		// force obsolete
		pose_left.pdb_info()->obsolete( true );

		// copy all residue information
		pose_left.pdb_info()->copy( *pose_right_.pdb_info(), 1, pose_right_.size(), original_left_nres + 1 );
	}

	assert( pose_left.size() == original_left_nres + pose_right_.size() );

	// END POSE MODIFY: after this point pose_left is static

	// store correct interval_
	interval_ = Interval( 1 + original_left_nres, pose_right_.size() + original_left_nres );

	// set the jump only if requested
	if ( use_rt_ ) {
		pose_left.conformation().set_stub_transform(
			StubID( core::pose::named_stub_id_to_stub_id( left_named_stub_id( left_position_ ) , pose_left ) ),
			StubID( core::pose::named_stub_id_to_stub_id( right_named_stub_id( right_position_ + original_left_nres  ), pose_left ) ),
			rt_
		);
	}
}


/// @brief do the actual reset of intervals, positions, etc to initial state
void ConnectRight::reset_accounting_impl() {
	left_position_ = original_interval().left;
	// right_position_ is never modified, so no reset necessary
	interval_ = Interval( 1, pose_right_.size() );
}


/// @brief initialization on construction
void ConnectRight::init() {
	left_stub_atoms_.push_back( "CA" );
	left_stub_atoms_.push_back( "N" );
	left_stub_atoms_.push_back( "CA" );
	left_stub_atoms_.push_back( "C" );

	right_stub_atoms_.push_back( "CA" );
	right_stub_atoms_.push_back( "N" );
	right_stub_atoms_.push_back( "CA" );
	right_stub_atoms_.push_back( "C" );

	assert( left_stub_atoms_.size() == 4 );
	assert( right_stub_atoms_.size() == 4 );
}


} // namespace build
} // namespace forge
} // namespace protocols

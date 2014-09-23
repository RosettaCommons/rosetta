// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/build/SegmentInsert.cc
/// @brief  insert an external segment flanked by new regions
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/SegmentInsert.hh>

// package headers
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>

// project headers
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/util.hh>
#include <core/kinematics/constants.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>


// numeric headers
#include <numeric/random/random.hh>

#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

#include <algorithm>


namespace protocols {
namespace forge {
namespace build {


// static for this file
static thread_local basic::Tracer TR( "protocols.forge.build.SegmentInsert" );


/// @brief default constructor
SegmentInsert::SegmentInsert() :
	Super()
{}


/// @brief sec.struct only constructor (poly-ala for flanking regions)
/// @param[in] interval The interval between which the insert will span.
///  To perform a pure insertion without replacing any residues
///  within a region, use an interval with a zero as the left endpoint, e.g.
///  [0, insert_after_this_residue].  If inserting before the first residue
///  of the Pose then interval = [0,0].  If inserting after the last residue
///  of the Pose then interval = [0, last_residue].
/// @param[in] ss The secondary structure specifying the flanking regions,
///  with a character '^' specifying where the insert is to be placed.
/// @param[in] insert The Pose to insert.
/// @param[in] keep_known_bb_torsions_at_junctions Attempt to keep the omega
///  at original_interval().left-1, the phi at original_interval().left, and
///  the psi+omega at original_interval().right present from the original Pose
///  in the modified Pose.  This should be false for pure insertions.
/// @param[in] connection_scheme Connect insertion on its N-side, C-side,
///  or decide randomly between the two (default RANDOM).
SegmentInsert::SegmentInsert(
	Interval const & i,
	String const & ss,
	Pose const & insert,
	bool const keep_known_bb_torsions_at_junctions,
	SegmentInsertConnectionScheme::Enum connection_scheme
) :
	Super( i, insert.residue( 1 ).residue_type_set().get_self_ptr() ),
	interval_( i ),
	ss_( ss ),
	keep_known_bb_torsions_at_junctions_( keep_known_bb_torsions_at_junctions ),
	insert_connection_scheme_( connection_scheme ),
	insert_pose_( insert )
{
	// construct the poly-ala a.a. string w/ insertion point
	for ( Size i = 0, ie = ss_.length(); i < ie; ++i ) {
		if ( ss.at( i ) != insertion_char() ) {
			aa_.push_back( 'A' );
		} else {
			aa_.push_back( insertion_char() );
		}
	}

	// safety
	if ( performing_pure_insertion() && keep_known_bb_torsions_at_junctions_ ) {
		TR.Warning << "keep_known_bb_torsions_at_junctions set to True, but performing pure insertion, so forcing the setting to False" << std::endl;
		keep_known_bb_torsions_at_junctions_ = false;
	}

	init();
}


/// @brief sec.struct + aa constructor
/// @param[in] interval The interval between which the insert will span.
///  To perform a pure insertion without replacing any residues
///  within a region, use an interval with a zero as the left endpoint, e.g.
///  [0, insert_after_this_residue].  If inserting before the first residue
///  of the Pose then interval = [0,0].  If inserting after the last residue
///  of the Pose then interval = [0, last_residue].
/// @param[in] ss The secondary structure specifying the flanking regions,
///  with a character '^' specifying where the insert is to be placed.
/// @param[in] aa The annotated amino acid string specifying the flanking
///  regions, with a character '^' specifying where the insert is to be
///  placed.
/// @param[in] insert The Pose to insert.
/// @param[in] keep_known_bb_torsions_at_junctions Attempt to keep the omega
///  at original_interval().left-1, the phi at original_interval().left, and
///  the psi+omega at original_interval().right present from the original Pose
///  in the modified Pose.  This should be false for pure insertions.
/// @param[in] connection_scheme Connect insertion on its N-side, C-side,
///  or decide randomly between the two (default RANDOM).
/// @remarks length of the *one-letter* aa must equal the length of ss
SegmentInsert::SegmentInsert(
	Interval const & i,
	String const & ss,
	String const & aa,
	Pose const & insert,
	bool const keep_known_bb_torsions_at_junctions,
	SegmentInsertConnectionScheme::Enum connection_scheme
) :
	Super( i, insert.residue( 1 ).residue_type_set().get_self_ptr() ),
	interval_( i ),
	ss_( ss ),
	aa_( aa ),
	keep_known_bb_torsions_at_junctions_( keep_known_bb_torsions_at_junctions ),
	insert_connection_scheme_( connection_scheme ),
	insert_pose_( insert )
{
	if ( !aa.empty() ) {
		// length and insertion point correspondence checks
		String const one_letter_aa = core::pose::annotated_to_oneletter_sequence( aa_ );
		runtime_assert( ss_.length() == one_letter_aa.length() );
		runtime_assert( ss_.find( insertion_char() ) == one_letter_aa.find( insertion_char() ) );

	} else {

		// construct the poly-ala a.a. string w/ insertion point
		for ( Size i = 0, ie = ss_.length(); i < ie; ++i ) {
			if ( ss.at( i ) != insertion_char() ) {
				aa_.push_back( 'A' );
			} else {
				aa_.push_back( insertion_char() );
			}
		}
	}

	// safety
	if ( performing_pure_insertion() && keep_known_bb_torsions_at_junctions_ ) {
		TR.Warning << "keep_known_bb_torsions_at_junctions set to True, but performing pure insertion, so forcing the setting to False" << std::endl;
		keep_known_bb_torsions_at_junctions_ = false;
	}

	init();
}


/// @brief copy constructor
SegmentInsert::SegmentInsert( SegmentInsert const & rval ) :
	Super( rval ),
	interval_( rval.interval_ ),
	ss_( rval.ss_ ),
	aa_( rval.aa_ ),
	insert_connection_scheme_( rval.insert_connection_scheme_ ),
	insert_pose_( rval.insert_pose_ ),
	insert_pose_torsion_override_movemap_( rval.insert_pose_torsion_override_movemap_ )
{}


/// @brief default destructor
SegmentInsert::~SegmentInsert() {}


/// @brief copy assignment
SegmentInsert & SegmentInsert::operator =( SegmentInsert const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		interval_ = rval.interval_;
		ss_ = rval.ss_;
		aa_ = rval.aa_;
		insert_connection_scheme_ = rval.insert_connection_scheme_;
		insert_pose_ = rval.insert_pose_;
		insert_pose_torsion_override_movemap_ = rval.insert_pose_torsion_override_movemap_;
	}

	return *this;
}


/// @brief clone this object
BuildInstructionOP SegmentInsert::clone() const {
	return BuildInstructionOP( new SegmentInsert( *this ) );
}


/// @brief the pose to insert
SegmentInsert::Pose const & SegmentInsert::insert_pose() const {
	return insert_pose_;
}


/// @brief get the absolute index of the insertion point with respect to the
///  flanking regions (i.e. the index inside the ss string)
/// @return the index, otherwise std::string::npos if not found
SegmentInsert::Size SegmentInsert::insertion_point_absolute_index() const {
	return ss_.find( insertion_char() );
}


/// @brief get the residue at the start of the insertion relative to the
///  modified interval (flanking positions are not part of the insertion!)
/// @return the residue position, otherwise 0 if not found
SegmentInsert::Size SegmentInsert::insertion_start_residue() const {
	return interval().left + ss_.find( insertion_char() );
}


/// @brief get the residue at the end of the insertion relative to the
///  modified interval (flanking positions are not part of the insertion!)
/// @return the residue position, otherwise 0 if not found
SegmentInsert::Size SegmentInsert::insertion_end_residue() const {
	return interval().left + ss_.find( insertion_char() ) + insert_pose_.n_residue() - 1;
}


/// @brief get the number of flanking residues to the left of the insertion
///  point
SegmentInsert::Size SegmentInsert::flanking_left_nres() const {
	return ss_.find( insertion_char() );
}


/// @brief get the number of flanking residues to the right of the insertion
///  point
SegmentInsert::Size SegmentInsert::flanking_right_nres() const {
	return ss_.length() - ss_.find( insertion_char() ) - 1;
}


/// @brief get the ss string of the flanking residues to the left of the
///  insertion point
SegmentInsert::String SegmentInsert::flanking_left_ss() const {
	return ss_.substr( 0, ss_.find( insertion_char() ) );
}


/// @brief get the ss string of the flanking residues to the right of the
///  insertion point
SegmentInsert::String SegmentInsert::flanking_right_ss() const {
	return ss_.substr( ss_.find( insertion_char() ) + 1 );
}


/// @brief get the annotated aa string of the flanking residues to the left
///  of the insertion point
SegmentInsert::String SegmentInsert::flanking_left_aa() const {
	return aa_.substr( 0, aa_.find( insertion_char() ) );
}


/// @brief get the annotated aa string of the flanking residues to the right
///  of the insertion point
SegmentInsert::String SegmentInsert::flanking_right_aa() const {
	return aa_.substr( aa_.find( insertion_char() ) + 1 );
}


/// @brief a copy of the working range of residues specifying the modified region
/// @details This residue range can change wrt length changes in Pose /Conformation
///  being watched.
Interval SegmentInsert::interval() const {
	return interval_;
}


/// @brief return a copy of the set of positions within the new region
///  that were pre-existing in the original Pose prior to modify()
/// @return An empty set, there are no pre-existing positions.
SegmentInsert::Positions SegmentInsert::preexisting_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions that are "new" and did
///  not exist in the original Pose.
SegmentInsert::Positions SegmentInsert::new_positions() const {
	using protocols::forge::methods::closed_range;

	// everything in the interval is new
	Interval const ival = interval();
	return closed_range( ival.left, ival.right );
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has a defined conformation.  E.g. existing or copied residues.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentInsert::Positions SegmentInsert::defined_positions() const {
	using protocols::forge::methods::closed_range;

	// only the insert is defined
	return closed_range( insertion_start_residue(), insertion_end_residue() );
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has an undefined conformation.  E.g. newly created residues.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentInsert::Positions SegmentInsert::undefined_positions() const {
	using protocols::forge::methods::insert_closed_range;

	Positions undefined;

	// the flanking regions are undefined
	Interval const ival = interval();
	insert_closed_range( ival.left, insertion_start_residue() - 1, undefined );
	insert_closed_range( insertion_end_residue() + 1, ival.right, undefined );

	return undefined;
}


/// @brief return a copy of the MoveMap that defines the moveable/fixed
///  positions/dofs for this instruction
/// @return TO BE FILLED IN
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentInsert::MoveMap SegmentInsert::movemap() const {
	using core::id::BB;
	using core::id::phi_torsion;
	using core::id::psi_torsion;
	using core::id::omega_torsion;
	using core::id::TorsionID;

	typedef MoveMap::MoveMapTorsionID MoveMapTorsionID;
	typedef MoveMap::TorsionTypeMap TorsionTypeMap;
	typedef MoveMap::MoveMapTorsionID_Map MoveMapTorsionID_Map;
	typedef MoveMap::TorsionID_Map TorsionID_Map;

	MoveMap mm;

	Interval const ival = interval();
	Interval const insertion( insertion_start_residue(), insertion_end_residue() );
	Interval const left_flank( ival.left, insertion.left - 1 );
	Interval const right_flank( insertion.right + 1, ival.right );

	// If attempting to keep the backbone torsions at in the positions equivalent
	// to the endpoints of the original interval in the original Pose, mark
	// partially moveable backbone at fixed-flanking and flanking-insertion
	// junction points
	if ( keep_known_bb_torsions_at_junctions_ ) {
		// 1. start of the left flanking region
		mm.set( TorsionID( left_flank.left, BB, phi_torsion ), false );
		mm.set( TorsionID( left_flank.left, BB, psi_torsion ), true );
		mm.set( TorsionID( left_flank.left, BB, omega_torsion ), true );

		// 2. left and right endpoints of the insertion
		mm.set( TorsionID( insertion.left, BB, phi_torsion ), true );
		mm.set( TorsionID( insertion.left, BB, psi_torsion ), false );
		mm.set( TorsionID( insertion.left, BB, omega_torsion ), false );
		mm.set( TorsionID( insertion.right, BB, phi_torsion), false );
		mm.set( TorsionID( insertion.right, BB, psi_torsion ), true );
		mm.set( TorsionID( insertion.right, BB, omega_torsion ), true );

		// 3. end of the right flanking region
		mm.set( TorsionID( right_flank.right, BB, phi_torsion ), true );
		mm.set( TorsionID( right_flank.right, BB, psi_torsion ), false );
		mm.set( TorsionID( right_flank.right, BB, omega_torsion ), false );

	} else { // all bb torsions moveable
		mm.set_bb( left_flank.left, true );
		mm.set_bb( insertion.left, true );
		mm.set_bb( insertion.left, true );
		mm.set_bb( right_flank.right, true );
	}

	// mark fully moveable parts of left flanking region
	for ( Size i = left_flank.left + 1; i <= left_flank.right; ++i ) {
		mm.set_bb( i, true );
	}

	// mark fully fixed parts of fixed insert
	for ( Size i = insertion.left + 1, ie = insertion.right - 1; i <= ie; ++i ) {
		mm.set_bb( i, false );
	}

	// mark fully moveable parts of right insert
	for ( Size i = right_flank.left, ie = right_flank.right - 1; i <= ie; ++i ) {
		mm.set_bb( i, true );
	}

	// Now import torsion settings in the insert_pose_torsion_override_movemap_,
	// transferring from lowest to highest stringency.
	Size const ridx_offset = insertion.left - 1; // residue index offset, add this to movemap settings to get proper numbering

	// TorsionType first
	for ( TorsionTypeMap::const_iterator i = insert_pose_torsion_override_movemap().torsion_type_begin(), ie = insert_pose_torsion_override_movemap().torsion_type_end(); i != ie; ++i ) {
		mm.set( i->first, i->second );
	}

	// MoveMapTorsionID second
	for ( MoveMapTorsionID_Map::const_iterator i = insert_pose_torsion_override_movemap().movemap_torsion_id_begin(), ie = insert_pose_torsion_override_movemap().movemap_torsion_id_end(); i != ie; ++i ) {
		MoveMapTorsionID shifted_id = i->first;
		shifted_id.first += ridx_offset;
		assert( insertion.left <= shifted_id.first && shifted_id.first <= insertion.right );
		mm.set( shifted_id, i->second );
	}

	// TorsionID last
	for ( TorsionID_Map::const_iterator i = insert_pose_torsion_override_movemap().torsion_id_begin(), ie = insert_pose_torsion_override_movemap().torsion_id_end(); i != ie ; ++i ) {
		TorsionID shifted_id = i->first;
		shifted_id.rsd() += ridx_offset;
		assert( insertion.left <= shifted_id.rsd() && shifted_id.rsd() <= insertion.right );
		mm.set( shifted_id, i->second );
	}

	return mm;
}


/// @brief set a torsion (bb/chi) specific override movemap indexed wrt the insert
///  Pose (residue indices may only be within the range [1, insert_pose.n_residue()]
/// @remarks When generating the movemap(), this torsion movemap will be enforced.
///  Only *explicit* settings of TorsionType, MoveMapTorsionID, and TorsionID will
///  be honored.  Implicit false settings are ignored.
void SegmentInsert::insert_pose_torsion_override_movemap( MoveMap const & mm ) {

	typedef MoveMap::TorsionTypeMap TorsionTypeMap;
	typedef MoveMap::MoveMapTorsionID_Map MoveMapTorsionID_Map;
	typedef MoveMap::TorsionID_Map TorsionID_Map;

	// run through all explicit torsion settings and make sure they're within
	// the range [1, insert_pose_.n_residue()] and only consist of BB & CHI
	// settings

	// TorsionTypeMap first
	for ( TorsionTypeMap::const_iterator i = mm.torsion_type_begin(), ie = mm.torsion_type_end(); i != ie; ++i ) {

		if ( i->first != core::id::BB && i->first != core::id::CHI ) {
			TR.Error << "ERROR: insert_pose_torsion_override_movemap() passed a MoveMap with an unhandled TorsionType setting" << std::endl;
			runtime_assert( false );
		}
	}

	// MoveMapTorsionID_Map second
	for ( MoveMapTorsionID_Map::const_iterator i = mm.movemap_torsion_id_begin(), ie = mm.movemap_torsion_id_end(); i != ie; ++i ) {

		if ( i->first.second != core::id::BB && i->first.second != core::id::CHI ) {
			TR.Error << "ERROR: insert_pose_torsion_override_movemap() passed a MoveMap with an unhandled"
				<< " MoveMapTorsionID type at residue " << i->first.first << std::endl;
			runtime_assert( false );
		}

		if ( i->first.first > insert_pose_.n_residue() ) {
			TR.Error << "ERROR: insert_pose_torsion_override_movemap() passed a MoveMap with a MoveMapTorsionID"
				<< " setting greater than the total number of residues in the insert pose at residue "
				<< i->first.first << std::endl;
			runtime_assert( false );
		}
	}

	// TorsionID_Map last
	for ( TorsionID_Map::const_iterator i = mm.torsion_id_begin(), ie = mm.torsion_id_end(); i != ie; ++i ) {

		if ( i->first.type() != core::id::BB && i->first.type() != core::id::CHI ) {
			TR.Error << "ERROR: insert_pose_torsion_override_movemap() passed a MoveMap with an unhandled"
				<< " TorsionID type at residue " << i->first.rsd() << std::endl;
			runtime_assert( false );
		}

		if ( i->first.rsd() > insert_pose_.n_residue() ) {
			TR.Error << "ERROR: insert_pose_torsion_override_movemap() passed a MoveMap with a MoveMapTorsionID"
				<< " setting greater than the total number of residues in the insert pose at residue "
				<< i->first.rsd() << std::endl;
			runtime_assert( false );
		}
	}

	// set the movemap
	insert_pose_torsion_override_movemap_ = mm;
}


/// @brief update indexing on residue append
void SegmentInsert::on_residue_append( LengthEvent const & event ) {
	if ( event.position < interval_.left ) {
		//++interval_.left;
		interval_.left += event.length_change;
	}

	if ( event.position < interval_.right ) {
		//++interval_.right;
		interval_.right += event.length_change;
	}
}


/// @brief update indexing on residue prepend
void SegmentInsert::on_residue_prepend( LengthEvent const & event ) {
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
void SegmentInsert::on_residue_delete( LengthEvent const & event ) {
	// event.position == interval.left is not caught below.
	// It has context dependent consequences and is manually corrected for
	// during modify().
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
SegmentInsert::Positions SegmentInsert::original_kept_positions() const {
	return Positions();
}


/// @brief return set of positions within the original interval that will
///  be deleted in this BuildInstruction
/// @return A set containing all positions in the open interval
///  (original.left, original.right)
SegmentInsert::Positions SegmentInsert::original_deleted_positions() const {
	using protocols::forge::methods::closed_range;
	return closed_range( original_interval().left, original_interval().right );
}


/// @brief return set of any fixed positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no fixed positions necessary
SegmentInsert::Positions SegmentInsert::original_fixed_positions() const {
	Positions fixed;

	fixed.insert( original_interval().left - 1 );
	fixed.insert( original_interval().right + 1 );

	return fixed;
}


/// @brief return set of any mutable positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
SegmentInsert::Positions SegmentInsert::original_mutable_positions() const {
	// for SegmentInsert, this is the same as the original deleted positions
	return original_deleted_positions();
}


/// @brief do the actual work of modifying the Pose
void SegmentInsert::modify_impl( Pose & pose ) {
	using core::chemical::ResidueTypeCOPs;
	using core::conformation::Residue;
	using core::conformation::ResidueOP;
	using core::conformation::ResidueOPs;
	using core::id::AtomID;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;
	using SegmentInsertConnectionScheme::C;
	using SegmentInsertConnectionScheme::N;
	using SegmentInsertConnectionScheme::RANDOM_SIDE;

	using core::conformation::get_anchor_atomno;
	using protocols::forge::methods::find_connecting_jump;
	using protocols::forge::methods::find_cutpoint;
	using protocols::forge::methods::grow_left_rtype;
	using protocols::forge::methods::grow_left_r;
	using protocols::forge::methods::grow_right_r;
	using protocols::forge::methods::grow_right_rtype;
	using protocols::forge::methods::order;
	using protocols::forge::methods::trans_omega;

	// there are two cases below:
	// (1) internal insertion where the removed segment may either be continuous
	//     or have cutpoints at any position along [original.left, original.right-1]
	// (2) boundary insertion where one of the endpoints of the removed segment
	//     is a terminus

	// are we doing an n- or c-terminal extension?
	bool performing_n_term_insertion = false;
	bool performing_c_term_insertion = false;

	if ( performing_pure_insertion() ) {
		performing_n_term_insertion = original_interval().right == 0;
		performing_c_term_insertion = original_interval().right == pose.n_residue();
	}
	assert( !( performing_n_term_insertion && performing_c_term_insertion ) );

	// Cache information from original structure.
	bool left_has_lower_terminus = false;
	bool left_has_upper_terminus = false;
	bool right_has_lower_terminus = false;
	bool right_has_upper_terminus = false;

	if ( !performing_n_term_insertion ) {

		Size const left_endpoint = performing_pure_insertion() ? interval_.right : interval_.left;
		left_has_lower_terminus = pose.residue( left_endpoint ).is_lower_terminus();
		left_has_upper_terminus = pose.residue( left_endpoint ).is_upper_terminus();

	} else { // case for fully n-terminal insertion
		left_has_lower_terminus = true;
		left_has_upper_terminus = false;
	}

	if ( !performing_c_term_insertion ) {

		Size const right_endpoint = performing_pure_insertion() ? interval_.right + 1 : interval_.right;
		right_has_lower_terminus = pose.residue( right_endpoint ).is_lower_terminus();
		right_has_upper_terminus = pose.residue( right_endpoint ).is_upper_terminus();

	} else { // case for fully c-terminal insertion
		right_has_lower_terminus = false;
		right_has_upper_terminus = true;
	}

	// If the following runtime_assert is triggered, that likely means you
	// tried to replace the entire length of the Pose or an entire chain of
	// a multi-chain Pose.  This class currently doesn't handle those cases
	// properly and will likely break.
	runtime_assert( !( left_has_lower_terminus && right_has_upper_terminus ) );

	// count # cutpoints along segment [left, right)
	Size n_cutpoints = 0;
	if ( performing_pure_insertion() ) {

		// for pure insertion/extensions we only check the residue
		// where the insertion will be made
		if ( pose.fold_tree().is_cutpoint( interval_.right ) ) {
			++n_cutpoints;
		}

	} else { // regular case

		for ( Size i = interval_.left; i < interval_.right; ++i ) {
			if ( pose.fold_tree().is_cutpoint( i ) ) {
				++n_cutpoints;
			}
		}
	}

	// save psi @ left-1 and phi @ right+1 to ensure they are set correctly
	// and not junked after residue deletion, set to 999.0 if not applicable
	Real const pre_psi = !( left_has_lower_terminus || left_has_upper_terminus ) && interval_.left > 1 ? pose.psi( interval_.left - 1 ) : 999.0;
	Real const post_phi = !( right_has_lower_terminus || right_has_upper_terminus ) && interval_.right < pose.n_residue() ? pose.phi( interval_.right + 1 ) : 999.0;

	// Store omega at original_interval().left-1, phi at original_interval().left,
	// and psi+omega at original_interval().right if user requests keeping them in
	// the equivalent positions in the modified Pose. Set to 999.0 if not applicable.
	//
	// TODO: I think there might be an issue here... if keeping the known
	// backbone torsions at the junctions, we should consider also keeping
	// in bond lengths and bond angles surrounding the junction.  That's
	// currently not done.
	Real left_endpoint_minus_one_omega = 999.0;
	Real left_endpoint_phi = 999.0;
	Real right_endpoint_psi = 999.0;
	Real right_endpoint_omega = 999.0;
	if ( keep_known_bb_torsions_at_junctions_ && !performing_pure_insertion() ) {
		// on the left
		if ( !left_has_lower_terminus && interval_.left > 1 ) {
			left_endpoint_minus_one_omega = pose.omega( interval_.left - 1 );
			left_endpoint_phi = pose.phi( interval_.left );
		}

		// on the right
		if ( !right_has_upper_terminus && interval_.right < pose.n_residue() ) {
			right_endpoint_psi = pose.psi( interval_.right );
			right_endpoint_omega = pose.omega( interval_.right );
		}
	}

	// grab residue types from aa string
	ResidueTypeCOPs r_types_flanking_left = core::pose::residue_types_from_sequence(
		flanking_left_aa(),
		residue_type_set(),
		false
	);

	ResidueTypeCOPs r_types_flanking_right = core::pose::residue_types_from_sequence(
		flanking_right_aa(),
		residue_type_set(),
		false
	);

	//testing
	//first identify the number of jups at this stage, this is the number to keep.
	Size num_jumps_pre_processing (pose.num_jump());

	// BEGIN INTERVAL SHIFT: after this point, interval_ will begin to shift due to length
	// changes in the Pose

	if ( !( performing_n_term_insertion || performing_c_term_insertion ) ) { // internal insertion

		// Two subcases:
		// (1) a pure insertion, no residue replacement
		// (2) insertion + residue replacement
		if ( !performing_pure_insertion() ) {
			// Blow away the old segment.  delete_residue_range_slow() should keep
			// all coordinates intact, though the bond lengths/angles adjacent to
			// the deletion might get funky.  Any strange geometry should get resolved
			// upon append/prepend of residues w/ ideal geometry below
			pose.conformation().delete_residue_range_slow( interval_.left, interval_.right );

			// At this point interval_.left should temporarily sit at left-1.
			--interval_.left;
		} else {
			interval_.left = interval_.right;
			++interval_.right;
		}

		// Step 0: create a cutpoint if necessary
		// Make a new jump only if the original segment was continuous.
		// If the original segment had jumps anywhere along [left, right)
		// this should imply there is already a jump connecting the left and right
		// halves.  If the jump existed somewhere in the original segment, the
		// residue deletion should have shifted the jump to an appropriate place.
		FoldTree ft = pose.fold_tree();
		Size const root = ft.root();

		if ( n_cutpoints == 0 && !left_has_lower_terminus && !right_has_upper_terminus
			&& !left_has_upper_terminus && !right_has_lower_terminus
		)
		{
			Size const new_jump_number = ft.new_jump( interval_.left, interval_.right, interval_.left );

			// following ensures changes in dihedral below do not shift the conformation
			ft.jump_edge( new_jump_number ).start_atom() = "N";
			ft.jump_edge( new_jump_number ).stop_atom() = "C";

		} else if ( n_cutpoints > 0 ) { // check/modify existing jump

			Size const jump = find_connecting_jump( interval_.left, interval_.right, ft );

			// We only do operations if there happens to be a jump defined
			// between the segments that 'left' and 'right' live on.
			if ( jump > 0 ) {

				Edge jump_edge = ft.jump_edge( jump );
				ft.delete_edge( jump_edge ); // remove edge

				// next two 'if' statements ensures changes in dihedral below do
				// not shift the conformation
				order( jump_edge );

				if ( static_cast< Size >( jump_edge.start() ) == interval_.left ) {
					jump_edge.start_atom() = "N";
				} else if ( jump_edge.start_atom().empty() ) {
					// if no existing jump atom, need to fill it in, otherwise
					// edge will not be recognized as having specific atom
					// settings
					jump_edge.start_atom() = pose.residue( jump_edge.start() ).atom_name(
						get_anchor_atomno( pose.residue( jump_edge.start() ), core::kinematics::dir_jump )
					);
				}

				if ( static_cast< Size >( jump_edge.stop() ) == interval_.right ) {
					jump_edge.stop_atom() = "C";
				} else if ( jump_edge.stop_atom().empty() ) {
					// if no existing jump atom, need to fill it in, otherwise
					// edge will not be recognized as having specific atom
					// settings
					jump_edge.stop_atom() = pose.residue( jump_edge.stop() ).atom_name(
						get_anchor_atomno( pose.residue( jump_edge.stop() ), core::kinematics::dir_jump )
					);
				}

				ft.add_edge( jump_edge ); // re-add modified edge
			}
		}

		ft.reorder( root );
		pose.fold_tree( ft );

	} else { // terminal insertion

		// need to manually modify the interval here for terminal insertions
		if ( performing_n_term_insertion ) {
			interval_.left = 1;
			interval_.right = 1;
		} else if ( performing_c_term_insertion ) {
			interval_.left = pose.n_residue();
			interval_.right = pose.n_residue();
		}

	}

	// Step 1: perform the growing of the flanking regions not directly
	// connected to the insertion
	Size before_insert_point = 0;
	if ( !left_has_lower_terminus ) { // grow towards the right
		assert( !performing_n_term_insertion );

		if ( !r_types_flanking_left.empty() ) {

			before_insert_point = grow_right_rtype(
				pose,
				interval_.left,
				r_types_flanking_left.begin(),
				r_types_flanking_left.end()
			);

		} else {
			before_insert_point = interval_.left;
		}
	}

	Size after_insert_point = 0;
	if ( !right_has_upper_terminus ) { // grow towards the left
		assert( !performing_c_term_insertion );

		if ( !r_types_flanking_right.empty() ) {

			after_insert_point = grow_left_rtype(
				pose,
				interval_.right,
				r_types_flanking_right.rbegin(),
				r_types_flanking_right.rend()
			);

		} else {
			after_insert_point = interval_.right;
		}
	}

	assert( !left_has_lower_terminus && !right_has_upper_terminus ? before_insert_point == after_insert_point - 1 : true );

	// Step 2: perform the insertion and perform any growing of any flanking
	// regions directly connected to the insertion.
	// The connection (peptide bond) between the junction of the insertion
	// and flanking regions will temporarily be distorted.  This is corrected
	// later in the procedure.

	// figure out exactly which side the insertion will be glued to
	SegmentInsertConnectionScheme::Enum insert_connection_scheme = insert_connection_scheme_;
	if ( left_has_lower_terminus ) {
		insert_connection_scheme = C;
	} else if ( right_has_upper_terminus ) {
		insert_connection_scheme = N;
	} else { // either choice is possible

		if ( insert_connection_scheme == RANDOM_SIDE ) {
			if ( numeric::random::rg().uniform() < 0.5 ) {
				insert_connection_scheme = N;
			} else {
				insert_connection_scheme = C;
			}
		}

	}

	assert( insert_connection_scheme != RANDOM_SIDE );

	// make a copy of the residues for safety
	ResidueOPs insert_residues;
	for ( Size i = 1, ie = insert_pose_.n_residue(); i <= ie; ++i ) {
		insert_residues.push_back( insert_pose_.residue( i ).clone() );
	}

	// add the residues from the insert Pose and grow any directly connected
	// flanking regions
	switch ( insert_connection_scheme ) {
		case N: {
			Size const new_anchor = grow_right_r(
				pose,
				before_insert_point, insert_residues.begin(), insert_residues.end(),
				true,
				true // place coordinates as-is
			);
			if ( right_has_upper_terminus ) {
				grow_right_rtype( pose, new_anchor, r_types_flanking_right.begin(), r_types_flanking_right.end() );
			}
			break;
		}
		case C: {
			Size const new_anchor = grow_left_r(
				pose,
				after_insert_point, insert_residues.rbegin(), insert_residues.rend(),
				true,
				true // place coordinates as-is
			);
			if ( left_has_lower_terminus ) {
				grow_left_rtype( pose, new_anchor, r_types_flanking_left.rbegin(), r_types_flanking_left.rend() );
			}
			break;
		}
		default:
			runtime_assert( false ); // should never get here
	}

	// finalize the interval so that it accurately reflects the new region
	if ( performing_n_term_insertion ) {
		interval_.left = 1;
		--interval_.right;
	} else if ( performing_c_term_insertion ) {
		++interval_.left;
		interval_.right = pose.n_residue();
	} else { // internal insertion {
		++interval_.left;
		--interval_.right;
	}

	assert( interval_.left < interval_.right ); // check ordering
	assert( interval_.length() > 0 );
	assert( interval_.length() == flanking_left_nres() + insert_pose_.n_residue() + flanking_right_nres() );

	// END INTERVAL SHIFT: after this point, interval_ has stabilized and stores
	// the new endpoints of the rebuilt segment

	// pick a cutpoint and set the new topology if doing internal insertion
	Size new_cutpoint = 0; // move this outside of the if statement, because the position is needed for editing no-jump tree
	if ( !( performing_n_term_insertion || performing_c_term_insertion ) ) {

		FoldTree new_ft = pose.fold_tree();
		Size const cutpoint = find_cutpoint( pose, interval_.left, interval_.right );
		assert( cutpoint > 0 ); // there should be a cutpoint

		switch ( insert_connection_scheme ) {
			case N: {
				// must remember to avoid making cut that affects bb torsions
				// at junction if user requests that option
				Size const largest_possible_cut_position = keep_known_bb_torsions_at_junctions_ ? interval_.right - 1 : interval_.right;

				if ( flanking_right_nres() > 0 ) {
					new_cutpoint = numeric::random::rg().random_range( interval_.right - flanking_right_nres(), largest_possible_cut_position );
				} else { // flanking_right_nres == 0
					new_cutpoint = interval_.right;
				}
				break;
			}
			case C: {
				if ( flanking_left_nres() > 0 ) {
					new_cutpoint = numeric::random::rg().random_range( interval_.left, interval_.left + flanking_left_nres() - 1 );
				} else {
					assert( interval_.left > 1 ); // safety, should never happen
					new_cutpoint = interval_.left - 1;
				}
				break;
			}
			default:
				runtime_assert( false ); // should never get here
		}

		if ( new_cutpoint != cutpoint ) {
			new_ft.slide_cutpoint( cutpoint, new_cutpoint );
		}

		// set the new topology
		pose.fold_tree( new_ft );
	}

	// Handle any jumps that happen to connect to left-1 or right+1.  These
	// need to have their jump atoms altered so that the torsion changes below
	// don't cause any rigid body shifts.
	FoldTree ft = pose.fold_tree();
	bool ft_was_altered = false;
	for ( FoldTree::const_iterator e = pose.fold_tree().begin(), ee = pose.fold_tree().end(); e != ee; ++e ) {
		if ( e->label() > 0 ) {
			Edge tmp = *e;
			order( tmp );

			bool start_changed = false;
			if ( static_cast< Size >( tmp.start() ) == interval_.left - 1 ) {
				tmp.start_atom() = "N";
				start_changed = true;
			}

			bool stop_changed = false;
			if ( static_cast< Size >( tmp.stop() ) == interval_.right + 1 ) {
				tmp.stop_atom() = "C";
				stop_changed = true;
			}

			// Must set opposing atom if it's not already set, otherwise
			// refold procedure will not recognize this edge as having
			// custom jump atoms.  Using default atom selection from
			// atom tree routines.
			if ( start_changed && !stop_changed && tmp.stop_atom().empty() ) {
				tmp.stop_atom() = pose.residue( tmp.stop() ).atom_name(
					get_anchor_atomno( pose.residue( tmp.stop() ), core::kinematics::dir_jump )
				);
				stop_changed = true;
			}

			if ( !start_changed && stop_changed && tmp.start_atom().empty() ) {
				tmp.start_atom() = pose.residue( tmp.start() ).atom_name(
					get_anchor_atomno( pose.residue( tmp.start() ), core::kinematics::dir_jump )
				);
				start_changed = true;
			}

			// re-add the edge only if the edge has been changed
			if ( start_changed && stop_changed ) {
				ft.delete_edge( *e );
				ft.add_edge( tmp );
				ft_was_altered = true;
			}
		}
	}

	if ( ft_was_altered ) {
		ft.reorder( pose.fold_tree().root() );
		pose.fold_tree( ft );
	}

 // special case for Remodel
  if (basic::options::option[basic::options::OptionKeys::remodel::no_jumps]){
    //ft.simple_tree(pose.total_residue());
    //idealize across the loop, in case they are not -- a big problem when taking out jumps
  //  protocols:loops::Loops loops_def_for_idealization;
   // loops_def_for_idealization.add_loop(protocols::loops::Loop(interval_.left-2,interval_.right+2, 1));

    utility::vector1<Size> cuts;
    utility::vector1< std::pair<Size,Size> > jumps;
    protocols::forge::methods::jumps_and_cuts_from_pose(pose,jumps, cuts);
    //debug
    //std::cout << "num cuts: " << cuts.size() << " num jumps: " << jumps.size() << std::endl;
    //translate:
    ObjexxFCL::FArray1D_int Fcuts( num_jumps_pre_processing);
    ObjexxFCL::FArray2D_int Fjumps(2, num_jumps_pre_processing);

  for ( Size i = 1; i<= num_jumps_pre_processing; ++i ) { // only keeping the old jumps, wipe new ones
      //std::cout << (int)jumps[i].first << " " << (int)jumps[i].second << std::endl;

      Fjumps(1,i) = std::min( (int)jumps[i].first, (int)jumps[i].second);
      Fjumps(2,i) = std::max( (int)jumps[i].first, (int)jumps[i].second);
      // DEBUG -- PRINT JUMPS AND CUTS
      //TR.Error << " jump " << i << " : " << jumps(1,i) << " , " << jumps(2,i) << std::endl;
    }
    //erase the new cut created
    cuts.erase( std::remove( cuts.begin(), cuts.end(), new_cutpoint ), cuts.end() );

    for ( Size i = 1; i<= cuts.size(); ++i ) {
      //std::cout << "cut " << (int)cuts[i] << std::endl;
      Fcuts(i) = (int)cuts[i];
  //    std::cout << " cut " << i << " : " << cuts(i) << std::endl;
    }

    // 4 make the foldtree
    FoldTree nojump_ft;
    nojump_ft.tree_from_jumps_and_cuts( pose.total_residue(), num_jumps_pre_processing, Fjumps, Fcuts, ft.root(), true ); // true );
    //std::cout << nojump_ft << std::endl;

    pose.fold_tree(nojump_ft);
    //std::cout << "IDEALIZE " << interval_.left << " to " << interval_.right << std::endl;
    //protocols::loops::set_extended_torsions_and_idealize_loops( pose, loops_def_for_idealization);
    for (Size i = interval_.left; i <= interval_.right; i++){
      core::conformation::idealize_position(i, pose.conformation());
    }
  }



	// Mark flanking regions and any existing connecting residues at left-1 and right
	// needing omega correction.  Omega correction is not actually performed until
	// bond length/angle corrections are done.
	Size omega_left_start = 0, omega_left_end = 0;
	Size omega_right_start = 0, omega_right_end = 0;

	if ( flanking_left_nres() > 0 ) {

		if ( interval_.left == 1 ) {
			omega_left_start = interval_.left;
			omega_left_end = flanking_left_nres();
		} else {
			omega_left_start = interval_.left - 1; // include residue immediately -1 to interval
			omega_left_end = omega_left_start + flanking_left_nres();
		}

		assert( omega_left_start <= omega_left_end );
	}

	if ( flanking_right_nres() > 0 ) {
		omega_right_start = interval_.right - flanking_right_nres();

		if ( interval_.right == pose.n_residue() ) {
			omega_right_end = interval_.right - 1;
		} else {
			omega_right_end = interval_.right;
		}

		assert( omega_right_start <= omega_right_end );
	}

	// Correct all bond angles and lengths at junction between flanking and
	// insert regions.  This will "rubber band snap" the insert into the right
	// position with existing (possibly non-ideal) internal geometry intact.
	switch ( insert_connection_scheme ) {
		case N:
			// always handle 'left_flanking <-> insert' junction residue, this
			// will never be a cutpoint
			if ( flanking_left_nres() > 0 ) {
				assert( omega_left_end > 0 );
				assert( !pose.fold_tree().is_cutpoint( omega_left_end ) );
				pose.conformation().insert_ideal_geometry_at_polymer_bond( omega_left_end );
			}

			// handle 'insert <-> right_flanking' junction residue only if not
			// a cutpoint
			if ( flanking_right_nres() > 0 ) {
				assert( omega_right_start > 0 );
				if ( !pose.fold_tree().is_cutpoint( omega_right_start ) ) {
					pose.conformation().insert_ideal_geometry_at_polymer_bond( omega_right_start );
				}
			}

			break;
		case C:
			// handle 'left_flanking <-> insert' junction residue only if not
			// a cutpoint
			if ( flanking_left_nres() > 0 ) {
				assert( omega_left_end > 0 );
				if ( !pose.fold_tree().is_cutpoint( omega_left_end ) ) {
					pose.conformation().insert_ideal_geometry_at_polymer_bond( omega_left_end );
				}
			}

			// always handle 'insert <-> right_flanking' junction residue, this
			// will never be a cutpoint
			if ( flanking_right_nres() > 0 ) {
				assert( omega_right_start > 0 );
				assert( !pose.fold_tree().is_cutpoint( omega_right_start ) );
				pose.conformation().insert_ideal_geometry_at_polymer_bond( omega_right_start );
			}

			break;
		default:
			runtime_assert( false ); // should never get here
	}

	// recover psi @ left-1 and phi @ right+1 to ensure they are set correctly
	// and not junked due to residue deletion, values set to 999.0 if not applicable
	if ( pre_psi <= 360.0 ) {
		pose.set_psi( interval_.left - 1, pre_psi );
		pose.set_omega( interval_.left - 1, 180.0 );
	}
	if ( post_phi <= 360.0 ) {
		pose.set_phi( interval_.right + 1, post_phi );
	}

	// correct all omegas for flanking regions and any existing connecting
	// residues at left-1 and right
	if ( flanking_left_nres() > 0 ) {
		trans_omega( omega_left_start, omega_left_end, pose );
	}

	if ( flanking_right_nres() > 0 ) {
		trans_omega( omega_right_start, omega_right_end, pose );
	}

	// recover old angles at junctions if requested
	if ( keep_known_bb_torsions_at_junctions_ && !performing_pure_insertion() ) {

		// torsion angles
		if ( left_endpoint_minus_one_omega <= 360.0 ) {
			pose.set_omega( interval_.left - 1, left_endpoint_minus_one_omega );
		}

		if ( left_endpoint_phi <= 360.0 ) {
			pose.set_phi( interval_.left, left_endpoint_phi );
		}

		if ( right_endpoint_psi <= 360.0 ) {
			pose.set_psi( interval_.right, right_endpoint_psi );
		}

		if ( right_endpoint_omega <= 360.0 ) {
			pose.set_omega( interval_.right, right_endpoint_omega );
		}

	}

	// set the desired secondary structure across the flanking regions and the
	// insert
	String const f_l_ss = flanking_left_ss();
	for ( Size i = 0, ie = f_l_ss.length(); i < ie; ++i ) {
		pose.set_secstruct( interval_.left + i, f_l_ss.at( i ) );
	}

	String const f_r_ss = flanking_right_ss();
	Size r = interval_.right - f_r_ss.length() + 1;
	for ( Size i = 0, ie = f_r_ss.length(); i < ie; ++i, ++r ) {
		pose.set_secstruct( r, f_r_ss.at( i ) );
	}

	r = interval_.left + flanking_left_nres();
	for ( Size i = 1, ie = insert_pose_.n_residue(); i <= ie; ++i, ++r ) {
		pose.set_secstruct( r, insert_pose_.secstruct( i ) );
	}

	// safety, make sure PDBInfo leaves obsolete
	if ( pose.pdb_info().get() ) {
		pose.pdb_info()->obsolete( true );
	}

}


/// @brief do the actual reset of intervals, positions, etc to initial state
void SegmentInsert::reset_accounting_impl() {
	interval_ = original_interval();
}


/// @brief init to be called during non-default constructors
void SegmentInsert::init() {
	using core::conformation::remove_lower_terminus_type_from_conformation_residue;
	using core::conformation::remove_upper_terminus_type_from_conformation_residue;

	// remove lower/upper terminus only at 1, nres
	if ( insert_pose_.n_residue() > 0 ) {
		if ( insert_pose_.residue( 1 ).is_lower_terminus() ) {
			core::conformation::remove_lower_terminus_type_from_conformation_residue( insert_pose_.conformation(), 1 );
		}

		if ( insert_pose_.residue( insert_pose_.n_residue() ).is_upper_terminus() ) {
			core::conformation::remove_upper_terminus_type_from_conformation_residue( insert_pose_.conformation(), insert_pose_.n_residue() );
		}
	}
}


} // namespace build
} // namespace forge
} // namespace protocols

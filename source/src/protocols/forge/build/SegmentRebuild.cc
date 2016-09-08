// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/SegmentRebuild.cc
/// @brief instruction to rebuild a segment
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/SegmentRebuild.hh>

// package headers
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>

// project headers

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/constants.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/PDBInfo.hh>


// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
// C++ headers
#include <algorithm>
#include <iostream>
#include <set>

#include <core/pose/annotated_sequence.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief default constructor
SegmentRebuild::SegmentRebuild() :
	Super(),
	keep_known_bb_torsions_at_junctions_( false )
{}


/// @brief sec.struct only constructor (poly-alanine for new region)
/// @param[in] interval rebuild this range of residues
/// @param[in] ss the secondary structure desired, also defines length of new build region
/// @param[in] rts the residue type set to use, default FA_STANDARD
/// @param[in] keep_known_bb_torsions_at_junctions Attempt to keep the omega
///  at original_interval().left-1, the phi at original_interval().left, and
///  the psi+omega at original_interval().right present from the original Pose
///  in the modified Pose.
/// @remarks length of the *one-letter* aa must equal the length of ss
SegmentRebuild::SegmentRebuild(
	Interval const & i,
	String const & ss,
	ResidueTypeSetCAP rts,
	bool const keep_known_bb_torsions_at_junctions
) :
	Super( i, rts ),
	interval_( i ),
	ss_( ss ),
	keep_known_bb_torsions_at_junctions_( keep_known_bb_torsions_at_junctions )
{
	// build poly-alanine if empty string
	if ( aa_.empty() ) {
		aa_ = String( ss_.length(), 'A' );
	}

	runtime_assert( ss_.length() == core::pose::annotated_to_oneletter_sequence( aa_ ).length() );
}


/// @brief full constructor
/// @param[in] interval rebuild this range of residues
/// @param[in] ss the secondary structure desired, also defines length of new build region
/// @param[in] aa the annotated amino acid sequence desired, default is poly-alanine
/// @param[in] rts the residue type set to use, default FA_STANDARD
/// @param[in] keep_known_bb_torsions_at_junctions Attempt to keep the omega
///  at original_interval().left-1, the phi at original_interval().left, and
///  the psi+omega at original_interval().right present from the original Pose
///  in the modified Pose.
/// @remarks length of the *one-letter* aa must equal the length of ss
SegmentRebuild::SegmentRebuild(
	Interval const & i,
	String const & ss,
	String const & aa,
	ResidueTypeSetCAP rts,
	bool const keep_known_bb_torsions_at_junctions
) :
	Super( i, rts ),
	interval_( i ),
	ss_( ss ),
	aa_( aa ),
	keep_known_bb_torsions_at_junctions_( keep_known_bb_torsions_at_junctions )
{
	// build poly-alanine if empty string
	if ( aa_.empty() ) {
		aa_ = String( ss_.length(), 'A' );
	}

	runtime_assert( ss_.length() == core::pose::annotated_to_oneletter_sequence( aa_ ).length() );
}


/// @brief copy constructor
SegmentRebuild::SegmentRebuild( SegmentRebuild const & rval ) :
	Super( rval ),
	interval_( rval.interval_ ),
	ss_( rval.ss_ ),
	aa_( rval.aa_ ),
	keep_known_bb_torsions_at_junctions_( rval.keep_known_bb_torsions_at_junctions_ )
{}


/// @brief default destructor
SegmentRebuild::~SegmentRebuild() {}


/// @brief copy assignment
SegmentRebuild & SegmentRebuild::operator =( SegmentRebuild const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		interval_ = rval.interval_;
		ss_ = rval.ss_;
		aa_ = rval.aa_;
		keep_known_bb_torsions_at_junctions_ = rval.keep_known_bb_torsions_at_junctions_;
	}
	return *this;
}


/// @brief clone this object
BuildInstructionOP SegmentRebuild::clone() const {
	return BuildInstructionOP( new SegmentRebuild( *this ) );
}


/// @brief return a copy of the set of positions within the new region
///  that were pre-existing in the original Pose prior to modify()
/// @return An empty set -- no positions are pre-existing.
SegmentRebuild::Positions SegmentRebuild::preexisting_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions that are "new" and did
///  not exist in the original Pose.
/// @return A set of positions spanning the entire modified interval -- all
///  positions are undefined.
SegmentRebuild::Positions SegmentRebuild::new_positions() const {
	using protocols::forge::methods::closed_range;

	Interval const ival = interval();
	return closed_range( ival.left, ival.right );
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has a defined conformation.  E.g. existing or copied residues.
/// @return An empty set -- no positions are defined.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentRebuild::Positions SegmentRebuild::defined_positions() const {
	return Positions();
}


/// @brief return a copy of the set of positions within the newly modified
///  region that has an undefined conformation.  E.g. newly created residues.
/// @return A set of positions spanning the entire modified interval -- all
///  positions are undefined.
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentRebuild::Positions SegmentRebuild::undefined_positions() const {
	// for SegmentRebuild this is the same as the new positions
	return new_positions();
}


/// @brief return a copy of the MoveMap that defines the moveable/fixed
///  positions/dofs for this instruction
/// @return a MoveMap with [interval.left, interval.right] bb & chi set to true
///  at the MoveMapTorsionID level
/// @details This set can change wrt length changes in Pose/Conformation being
///  watched.
SegmentRebuild::MoveMap SegmentRebuild::movemap() const {
	using core::id::BB;
	using core::id::phi_torsion;
	using core::id::psi_torsion;
	using core::id::omega_torsion;
	using core::id::TorsionID;

	MoveMap mm;

	Interval const ival = interval();
	Size begin = ival.left;
	Size end = ival.right;

	if ( keep_known_bb_torsions_at_junctions_ ) {
		++begin;
		--end;

		// left residue, only psi+omega+chi move
		mm.set( TorsionID( ival.left, BB, phi_torsion ), false );
		mm.set( TorsionID( ival.left, BB, psi_torsion ), true );
		mm.set( TorsionID( ival.left, BB, omega_torsion ), true );
		mm.set_chi( ival.left, true );

		// right residue, only phi+chi moves
		mm.set( TorsionID( ival.right, BB, phi_torsion ), true );
		mm.set( TorsionID( ival.right, BB, psi_torsion ), false );
		mm.set( TorsionID( ival.right, BB, omega_torsion ), false );
		mm.set_chi( ival.right, true );
	}

	for ( Size i = begin; i <= end; ++i ) {
		mm.set_bb( i, true );
		mm.set_chi( i, true );
	}

	return mm;
}


/// @brief update indexing on residue append
/// @remarks left and right endpoints of the interval can travel independently
void SegmentRebuild::on_residue_append( LengthEvent const & event ) {
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
/// @remarks left and right endpoints of the interval can travel independently
void SegmentRebuild::on_residue_prepend( LengthEvent const & event ) {
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
/// @remarks Left and right endpoints of the interval can travel independently.
void SegmentRebuild::on_residue_delete( LengthEvent const & event ) {
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
/// @return An empty set -- no positions are kept.
SegmentRebuild::Positions SegmentRebuild::original_kept_positions() const {
	return Positions();
}


/// @brief return set of positions within the original interval that will
///  be deleted in this BuildInstruction
/// @return A set containing all positions in the original interval.
SegmentRebuild::Positions SegmentRebuild::original_deleted_positions() const {
	using protocols::forge::methods::closed_range;
	return closed_range( original_interval().left, original_interval().right );
}


/// @brief return set of any fixed positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.  If SegmentRebuild has dependencies,
///  the set of positions shrinks from the endpoints of [original_left - 1, original_right + 1]
///  down to an empty set in the assumption that the user is requesting an
///  advanced feature and knows what they're doing around the endpoints,
///  e.g. rebuilding directly adjacent to a swapped out section.
/// @return empty set if no fixed positions necessary
SegmentRebuild::Positions SegmentRebuild::original_fixed_positions() const {
	Positions fixed;

	// If there are dependencies, then we relax the requirements and assume
	// that the user might be directing the rebuild section directly adjacent to
	// an e.g. swapped section and return an empty set.  This is an advanced
	// option, so we assume the user understands what they're doing around the
	// endpoints.

	if ( !has_dependencies() ) { // typical case without dependencies
		fixed.insert( original_interval().left - 1 );
		fixed.insert( original_interval().right + 1 );
	}

	return fixed;
}


/// @brief return set of any mutable positions necessary with respect to the original
///  interval and original Pose numbering
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.
/// @return empty set if no mutable positions
SegmentRebuild::Positions SegmentRebuild::original_mutable_positions() const {
	using protocols::forge::methods::closed_range;
	return closed_range( original_interval().left, original_interval().right );
}


/// @brief do the actual work of modifying the Pose
/// @details There are two cases:
///  <ul>
///        <li> internal rebuild where no residue of the segment is a terminus
///        <li> boundary rebuild where one residue of the segment is a terminus
///  </ul>
///  In the case of an internal rebuild, the old segment will be removed and
///  a new segment containing a random cutpoint will be inserted.  The new jump
///  will be created between (left-1) and (right+1) of the new segment.
void SegmentRebuild::modify_impl( Pose & pose ) {
	using core::chemical::ResidueTypeCOPs;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	using core::conformation::get_anchor_atomno;
	using protocols::forge::methods::closest_larger_peptide_vertex;
	using protocols::forge::methods::closest_smaller_peptide_vertex;
	using protocols::forge::methods::find_connecting_jump;
	using protocols::forge::methods::grow_left_rtype;
	using protocols::forge::methods::grow_right_rtype;
	using protocols::forge::methods::jumps_connected_to_position;
	using protocols::forge::methods::remove_cutpoint_variants;
	using protocols::forge::methods::order;
	using protocols::forge::methods::trans_omega;
	using namespace core::conformation;
	/*

	//special case for two chain build -- danger, primitive at the moment, only
	//works with non-N,C term extension.  Try placing in SegRebuld
	if (basic::options::option[basic::options::OptionKeys::remodel::two_chain_tree]()){
	Size second_start = basic::options::option[basic::options::OptionKeys::remodel::two_chain_tree];
	Size nres( pose.size());

	//FoldTree f(pose.fold_tree());
	FoldTree f;
	//make cutpoint
	f.add_edge(1, second_start-1, Edge::PEPTIDE);
	f.add_edge(second_start, nres, Edge::PEPTIDE);
	f.add_edge(second_start-1,second_start,1);//jump across the cut
	f.reorder(nres);
	pose.fold_tree(f);

	protocols:loops::Loops chain_def_loops;
	chain_def_loops.add_loop(protocols::loops::Loop(1,second_start));
	chain_def_loops.add_loop(protocols::loops::Loop(second_start,pose.size()));

	//add virtual residue, as star foldtree requires it
	if (pose.residue( pose.size()).aa() != core::chemical::aa_vrt){
	pose.append_residue_by_jump(*core::conformation::ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map("VRT")), pose.size());
	}

	//update foldtree to new foldtree
	f = pose.fold_tree();

	//update nres to include the new residue
	nres =  pose.size();

	f.reorder(nres);
	pose.fold_tree(f);
	protocols::forge::methods::make_star_foldtree(pose, chain_def_loops);
	}
	*/
	//typedef utility::vector1< Edge > Edges;

	//for nojump operation, need to keep track of the new cuts introduced
	Size new_cut = 0 ;

	// there are two cases below:
	// (1) internal rebuild where the segment may either be continuous or have cutpoints
	//     at any position along [left-1, right]
	// (2) boundary rebuild where one of the endpoints of the segment is a terminus

	// cache information from original structure
	bool const has_lower_terminus = pose.residue( interval_.left ).is_lower_terminus();
	bool const has_upper_terminus = pose.residue( interval_.right ).is_upper_terminus();
	Size const old_root = pose.fold_tree().root();

	// If the following runtime_assert is triggered, that likely means you
	// tried to replace the entire length of the Pose or an entire chain of
	// a multi-chain Pose.  This class currently doesn't handle those cases
	// properly and will likely break.
	// but when only the case of single chain is allowed
	if ( pose.conformation().num_chains() != 1 ) {
		runtime_assert( !( has_lower_terminus && has_upper_terminus ) );
	}

	// count # cutpoints along segment [left-1, right]
	Size n_cutpoints = 0;
	for ( Size i = interval_.left - 1; i <= interval_.right; ++i ) {
		if ( pose.fold_tree().is_cutpoint( i ) ) {
			++n_cutpoints;
		}
	}

	// save psi @ left-1 and phi @ right+1 to ensure they are set correctly
	// and not junked after residue deletion, set to 999.0 if not applicable
	Real const pre_psi = !has_lower_terminus ? pose.psi( interval_.left - 1 ) : 999.0;
	Real const post_phi = !has_upper_terminus ? pose.phi( interval_.right + 1 ) : 999.0;

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
	if ( keep_known_bb_torsions_at_junctions_ ) {
		// on the left
		if ( !has_lower_terminus && interval_.left > 1 ) {
			left_endpoint_minus_one_omega = pose.omega( interval_.left - 1 );
			left_endpoint_phi = pose.phi( interval_.left );
		}

		// on the right
		if ( !has_upper_terminus && interval_.right < pose.size() ) {
			right_endpoint_psi = pose.psi( interval_.right );
			right_endpoint_omega = pose.omega( interval_.right );
		}
	}

	// grab residue types from aa string
	ResidueTypeCOPs r_types = core::pose::residue_types_from_sequence( aa_, residue_type_set(), false );
	assert( ss_.length() == r_types.size() );

	// BEGIN INTERVAL SHIFT: after this point, interval_ will begin to shift due to length
	// changes in the Pose

	if ( pose.conformation().num_chains() == 1 && has_lower_terminus && has_upper_terminus ) {

		// Store observer state.  If this object was attaached to the Pose, we need
		// to re-attach.
		bool const was_connected = length_obs_link().valid();

		// Catch special case for single chain.  FoldTree currently refuses to delete
		// an entire tree, so we clear the entire Pose.
		pose.clear();

		// Clearing the Pose invalidates observers, so we need to re-hook ourselves
		// back onto the Pose if previously connected.
		if ( was_connected ) {
			attach_to( pose );
		}

	} else {

		// Blow away the old segment.  delete_residue_range_slow() should keep
		// all coordinates intact, though the bond lengths/angles adjacent to
		// the deletion might get funky.  Any strange geometry should get resolved
		// upon append/prepend of residues w/ ideal geometry below
		pose.conformation().delete_residue_range_slow( interval_.left, interval_.right );
		assert( interval_.left == interval_.right );

	}

	//testing
	//first identify the number of jups at this stage, this is the number to keep.
	Size num_jumps_pre_processing (pose.num_jump());


	// perform the append/prepend operations
	if ( pose.conformation().num_chains() == 1 && has_lower_terminus && has_upper_terminus ) {

		--interval_.left;
		interval_.right = 1; // artificially set to 1 for assertion check

		// only grow right; termini will not be recovered here due to prior
		// residue deletions and is corrected for below
#ifndef NDEBUG
		Size const right_endpoint = grow_right_rtype( pose, interval_.left, r_types.begin(), r_types.end() );
#else
		grow_right_rtype( pose, interval_.left, r_types.begin(), r_types.end() );
#endif
		assert( right_endpoint == interval_.right - 1 );

		// set to proper interval
		assert( ss_.length() == interval_.right - 1);
		interval_.left = 1;
		interval_.right = ss_.length(); // or --interval_.right

		// re-add terminus
		if ( !pose.residue( interval_.right ).is_upper_terminus() ) {
			core::pose::add_upper_terminus_type_to_pose_residue( pose, interval_.right );
		}

		// re-add terminus
		if ( !pose.residue( interval_.left ).is_lower_terminus() ) {
			core::pose::add_lower_terminus_type_to_pose_residue( pose, interval_.left );
		}

	} else if ( has_lower_terminus ) { // left boundary rebuild

		// At this point the endpoints of the interval should both sit at
		// the anchor position where the prepend operations should happen.

		// only grow left; termini will not be recovered here due to prior
		// residue deletions and is corrected for below
		Size const left_endpoint = grow_left_rtype( pose, interval_.right, r_types.rbegin(), r_types.rend() );
		assert( left_endpoint == interval_.right - r_types.size() );

		// correct endpoints of interval to match the actual rebuilt segment
		interval_.left = left_endpoint;
		--interval_.right;
		assert( interval_.right - interval_.left + 1 == r_types.size() );

		// re-add terminus
		if ( !pose.residue( interval_.left ).is_lower_terminus() ) {
			core::pose::add_lower_terminus_type_to_pose_residue( pose, interval_.left );
		}

		// re-root tree towards the right
		if ( original_interval().contains( old_root ) ) {
			FoldTree ft = pose.fold_tree();
			ft.reorder( closest_larger_peptide_vertex( interval_.right, ft, 2 ) );
			pose.fold_tree( ft );
		}

	} else if ( has_upper_terminus ) { // right boundary rebuild

		// At this point the endpoints of the interval should both sit +1
		// position to where the append operations should happen.
		// Temporarily shift both endpoints back one to sit at anchor.
		--interval_.left;
		--interval_.right;

		// only grow right; termini will not be recovered here due to prior
		// residue deletions and is corrected for below
		Size const right_endpoint = grow_right_rtype( pose, interval_.left, r_types.begin(), r_types.end() );
		assert( right_endpoint == interval_.left + r_types.size() );

		// correct endpoints of interval to match the actual rebuilt segment
		++interval_.left;
		interval_.right = right_endpoint;
		assert( interval_.right - interval_.left + 1 == r_types.size() );

		// re-add terminus
		if ( !pose.residue( interval_.right ).is_upper_terminus() ) {
			core::pose::add_upper_terminus_type_to_pose_residue( pose, interval_.right );
		}

		// re-root tree towards the left
		if ( original_interval().contains( old_root ) ) {
			FoldTree ft = pose.fold_tree();
			ft.reorder( closest_smaller_peptide_vertex( interval_.left, ft, 2 ) );
			pose.fold_tree( ft );
		}

	} else { // internal rebuild

		// assign cutpoint, -1 ensures the cut happens internal to the range
		// TODO: consider using sec.struct to preferentially select a cutpoint
		Size cut_index = numeric::random::rg().random_range( 1, r_types.size() - 1);

		//in order to synchronize all cutting, unfortunately this has to be in this
		//file
		if ( basic::options::option[basic::options::OptionKeys::remodel::RemodelLoopMover::bypass_closure]() ) {
			if ( basic::options::option[basic::options::OptionKeys::remodel::RemodelLoopMover::force_cutting_N]() ) {
				cut_index = 1;
			} else if ( basic::options::option[basic::options::OptionKeys::remodel::RemodelLoopMover::force_cutting_index].user() ) {
				cut_index = basic::options::option[basic::options::OptionKeys::remodel::RemodelLoopMover::force_cutting_index];
			} else {
				cut_index = r_types.size()-1;
			}
		}

		// Temporarily shift the endpoints to the anchor points:
		// The left endpoint goes back one and the right endpoint is
		// left where it is, so it can travel as the append/prepend operations
		// occur.
		--interval_.left;
		assert( interval_.left + 1 == interval_.right );

		// Do any necessary chain corrections before growing residues; need to
		// merge chains if different.  Here we check the position currently at
		// interval.left - 1 and remove it from the chain endings list.
		if (
				std::find(
				pose.conformation().chain_endings().begin(),
				pose.conformation().chain_endings().end(),
				interval_.left
				) != pose.conformation().chain_endings().end()
				) {
			pose.conformation().delete_chain_ending( interval_.left );
		}

		// Make a new jump only if the original segment was continuous.
		// If the original segment had jumps anywhere along [left-1, right]
		// this should imply there is already a jump connecting the left and right
		// halves.  If the jump existed somewhere in the original segment, the
		// residue deletion should have shifted the jump to an appropriate place.
		// In the case of an existing jump that is adjacent to the interval,
		// we need to modify the jump atoms so that changes in the dihedrals do
		// not shift the conformation.
		FoldTree ft = pose.fold_tree();
		Size const root = ft.root();

		if ( n_cutpoints == 0 ) { // need a new jump
			Size const new_jump_number = ft.new_jump( interval_.left, interval_.right, interval_.left );

			// following ensures changes in dihedral below do not shift the
			// conformation
			ft.jump_edge( new_jump_number ).start_atom() = "N";
			ft.jump_edge( new_jump_number ).stop_atom() = "C";

		} else { // check/modify existing jump

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

		// Append up until & including cutpoint.  Since this is an internal
		// rebuild, tell grow_right_rtype not to maintain any existing termini
		// variant; this is a safety, since residue deletion should have removed
		// any such variant.
		Size left_of_cut = grow_right_rtype(
			pose, interval_.left,
			r_types.begin(), r_types.begin() + cut_index,
			false
		);
		new_cut = left_of_cut;
		//std::cout << "new cuts made: " << left_of_cut << std::endl;
		// Prepend down until right before cutpoint.  Since this is an internal
		// rebuild, tell grow_left_rtype not to maintain any existing termini
		// variant; this is a safety, since residue deletion should have removed
		// any such variant.
		Size right_of_cut = grow_left_rtype(
			pose, interval_.right,
			r_types.rbegin(), r_types.rbegin() + ( r_types.size() - cut_index ),
			false
		);

		runtime_assert( left_of_cut + 1 == right_of_cut );
		runtime_assert( !pose.residue( left_of_cut ).is_upper_terminus() );
		runtime_assert( !pose.residue( right_of_cut ).is_lower_terminus() );

		// correct endpoints of interval to match the actual rebuilt segment
		++interval_.left;
		--interval_.right;
		assert( interval_.right - interval_.left + 1 == r_types.size() );
	}

	// END INTERVAL SHIFT: after this point, interval_ has stabilized and stores
	// the new endpoints of the rebuilt segment

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
	//std::cout << ft << std::endl;

	// special case for Remodel
	if ( basic::options::option[basic::options::OptionKeys::remodel::no_jumps] ) {
		//ft.simple_tree(pose.size());
		//idealize across the loop, in case they are not -- a big problem when taking out jumps
		// protocols:loops::Loops loops_def_for_idealization;
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
		cuts.erase( std::remove( cuts.begin(), cuts.end(), new_cut ), cuts.end() );

		for ( Size i = 1; i<= cuts.size(); ++i ) {
			//std::cout << "cut " << (int)cuts[i] << std::endl;
			Fcuts(i) = (int)cuts[i];
			//  std::cout << " cut " << i << " : " << cuts(i) << std::endl;
		}

		// 4 make the foldtree
		FoldTree nojump_ft;

		if ( basic::options::option[basic::options::OptionKeys::remodel::reroot_tree].user() ) { //rerooting tree so can anchor in stationary part of structure.
			Size new_root =  basic::options::option[basic::options::OptionKeys::remodel::reroot_tree];
			nojump_ft.tree_from_jumps_and_cuts( pose.size(), num_jumps_pre_processing, Fjumps, Fcuts, new_root , true ); // true );
		} else { //default
			nojump_ft.tree_from_jumps_and_cuts( pose.size(), num_jumps_pre_processing, Fjumps, Fcuts, ft.root(), true ); // true );
		}

		//std::cout << nojump_ft << std::endl;

		pose.fold_tree(nojump_ft);
		//std::cout << "IDEALIZE " << interval_.left << " to " << interval_.right << std::endl;
		//protocols::loops::set_extended_torsions_and_idealize_loops( pose, loops_def_for_idealization);
		for ( Size i = interval_.left; i <= interval_.right; i++ ) {
			core::conformation::idealize_position(i, pose.conformation());
		}
		//pose.dump_pdb("test_idl.pdb");
	}


	// recover psi @ left-1 and phi @ right+1 to ensure they are set correctly
	// and not junked due to residue deletion, values set to 999.0 if not applicable
	if ( pre_psi <= 360.0 ) {
		pose.set_psi( interval_.left - 1, pre_psi );
	}
	if ( post_phi <= 360.0 ) {
		pose.set_phi( interval_.right + 1, post_phi );
	}

	// assume proper omega
	Size const omega_left = ( interval_.left == 1 ) ? interval_.left : interval_.left - 1;
	Size const omega_right = ( interval_.right == pose.size() ) ? interval_.right - 1 : interval_.right;
	trans_omega( omega_left, omega_right, pose );

	// recover old angles at junctions if requested
	if ( keep_known_bb_torsions_at_junctions_ ) {

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

	// set the desired secondary structure
	for ( Size r = interval_.left, i = 0; r <= interval_.right; ++r, ++i ) {
		pose.set_secstruct( r, ss_.at( i ) );
	}

	// safety, make sure PDBInfo leaves obsolete
	if ( pose.pdb_info().get() ) {
		pose.pdb_info()->obsolete( true );
	}


}


/// @brief do the actual reset of intervals, positions, etc to initial state
void SegmentRebuild::reset_accounting_impl() {
	interval_ = original_interval();
}


} // namespace build
} // namespace forge
} // namespace protocols

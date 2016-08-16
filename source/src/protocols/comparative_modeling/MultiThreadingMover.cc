// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/MultiThreadingMover.hh
/// @brief
/// @author James Thompson

#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/comparative_modeling/MultiThreadingMover.hh>

#include <core/types.hh>


#include <core/chemical/ResidueTypeSet.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/Residue.fwd.hh>

#include <core/fragment/FragSet.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>


#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <basic/Tracer.hh>


#include <utility/vector1.hh>

// C++ headers
#include <string>

//Auto Headers
#include <core/pose/util.tmpl.hh>
// option key includes


namespace protocols {
namespace comparative_modeling {

MultiThreadingMover::MultiThreadingMover(
	Alignments const & alignments,
	Poses const & template_poses
) :
	protocols::moves::Mover( "MultiThreadingMover" ),
	template_poses_( template_poses ),
	alignments_( alignments ),
	build_query_loops_( true ),
	repack_query_( true ),
	randomize_loop_coords_( false ),
	min_loop_size_( 3 )
{}

//boolean setters
void MultiThreadingMover::build_loops( bool setting ) {
	build_query_loops_ = setting;
}

void MultiThreadingMover::randomize_loop_coords( bool setting ) {
	randomize_loop_coords_ = setting;
}

void MultiThreadingMover::repack_query( bool setting ) {
	repack_query_ = setting;
}

//boolean getters
bool MultiThreadingMover::build_loops() const {
	return build_query_loops_;
}

bool MultiThreadingMover::repack_query() const {
	return repack_query_;
}

bool MultiThreadingMover::randomize_loop_coords() {
	return randomize_loop_coords_;
}

void MultiThreadingMover::min_loop_size( core::Size const new_size ) {
	min_loop_size_ = new_size;
}

core::Size MultiThreadingMover::min_loop_size() const {
	return min_loop_size_;
}

utility::vector1< core::fragment::FragSetOP > MultiThreadingMover::frag_libs() const {
	return frag_libs_;
}

void MultiThreadingMover::frag_libs(
	utility::vector1< core::fragment::FragSetOP > new_libs
) {
	frag_libs_ = new_libs;
}

void MultiThreadingMover::apply(
	core::pose::Pose & query_pose
) {
	using std::set;
	using core::Size;
	using core::pose::Pose;
	using core::sequence::SequenceAlignment;
	static basic::Tracer tr(
		"protocols.comparative_modeling.MultiThreadingMover.apply"
	);

	check_internals();

	// loop over alignments and template poses, thread from each
	// one.

	core::id::AtomID_Mask missing( true );
	core::pose::initialize_atomid_map( missing, query_pose ); // used for repacking atoms
	for ( Size idx = 1; idx <= template_poses_.size(); ++idx ) {
		Pose const & template_pose( template_poses_[idx] );
		SequenceAlignment const & aln( alignments_[idx] );

		ThreadingMover threader( aln, template_pose );
		threader.build_loops(false);
		threader.repack_query(false);
		threader.randomize_loop_coords(false);
		threader.apply(query_pose);
		core::id::SequenceMapping map( aln.sequence_mapping(1,2) );
		for ( Size resi = 1; resi <= query_pose.total_residue(); ++resi ) {
			if ( map[resi] != 0 ) {
				for ( Size atomj = 1; atomj <= query_pose.residue(resi).natoms(); ++atomj ) {
					missing[ core::id::AtomID( atomj, resi ) ] = false;
				}
			}
		}
	}
	randomize_selected_atoms( query_pose, missing );

	// rebuild loops

	// repack query
} // apply

std::string
MultiThreadingMover::get_name() const {
	return "MultiThreadingMover";
}

void MultiThreadingMover::check_internals() const {
	runtime_assert( template_poses_.size() == alignments_.size() );
}

} // comparative_modeling
} // protocols

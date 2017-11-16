// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <devel/cycles/ReindexingMover.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>

namespace devel {
namespace cycles {

using namespace core;

static basic::Tracer tr( "devel.cycles.ReindexingMover" );

ReindexingMover::ReindexingMover(Size offset) { offset_ = offset; }

ReindexingMover::~ReindexingMover() = default;

/// The offset parameter passed to the constructor is used to determine the
/// displacement between the input and output poses.  The given pose should
/// have been passed through the SetupMover at some point, otherwise this
/// method may trigger some assertion failures.

void ReindexingMover::apply(pose::Pose &source_pose) {

	// Both the index and residue variable must be signed integers, because the
	// offset can be negative.  The remainder operator (%) behaves incorrectly
	// when one operand is signed and the other is unsigned.  This is noteworthy
	// because it is more common within the rosetta codebase to use the Size
	// type, which is unsigned.

	pose::Pose target_pose;
	int residues = source_pose.size();

	for ( int i = 0; i < residues; i++ ) {

		// There are a few things to note about how the offset index is calculated.
		// First, the for-loop index variable counts from 0 in order to facilitate
		// the remainder arithmetic.  The index variable is subsequently modified
		// to count from 1, in order to agree with the rest of rosetta.  If the
		// resulting index is less than one, it is modified again to move it back
		// into the range (1, residues).

		int index = 1 + (i - offset_) % residues;
		if ( index < 1 ) index += residues;

		conformation::Residue residue = source_pose.residue(index);
		target_pose.append_residue_by_bond(residue);
	}

	// It is necessary to clear to source pose when it has the same sequence as
	// the target pose and offset != 0.  This is because the assignment operator
	// in the pose class takes some shortcuts when it notices that the target has
	// the same length and sequence as the source.  In this case, these shortcuts
	// are illegal and end up crashing the program.  Clearing the pose side-steps
	// this problem.

	source_pose.clear();
	source_pose = target_pose;
}

} // End 'cycles' namespace.
} // End 'devel' namespace.


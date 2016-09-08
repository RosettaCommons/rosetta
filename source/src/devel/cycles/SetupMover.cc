// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <devel/cycles/SetupMover.hh>
#include <devel/cycles/ReindexingMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/NamedAtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>

namespace devel {
namespace cycles {

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "devel.cycles.SetupMover" );

SetupMover::SetupMover() {}

SetupMover::~SetupMover() = default;

std::string SetupMover::get_name() const {
	return "SetupMover"; 
}

void SetupMover::apply(pose::Pose& pose) {
	remove_termini_patches(pose);
	restore_amide_hydrogen(pose);
}

void SetupMover::remove_termini_patches(pose::Pose& pose) {

	Size first_index = 1;
	Size last_index = pose.size();

	// This method triggers three error messages in the Residue class.  Since all 
	// of these errors are corrected later on by this mover, these messages are 
	// needlessly worrisome.  There's no way to specifically mute the Residue 
	// tracer, so it's necessary to temporarily mute all tracers globally.
	//
	// Note that the import_pose::pose_from_file() method also generates an error 
	// message when loading cyclic peptides.  Unfortunately, there's nothing that 
	// can be done about that from within this mover.

	TR.super_mute(true);

	pose::remove_variant_type_from_pose_residue(
			pose, "LOWER_TERMINUS", first_index);

	TR.super_mute(false);

	pose::remove_variant_type_from_pose_residue(
			pose, "UPPER_TERMINUS", last_index);

}

/// When the N-terminal patch is removed, rosetta doesn't know how to rebuild 
/// the amide hydrogen because there's no residue on the far side of that amide 
/// bond.  This hydrogen can be rebuilt once this residue exists, which happens 
/// when the cycle is reindexed.  This method rotates the cycle forward on 
/// residue, rebuilds the missing hydrogen, then rotates the cycle back to 
/// where it started.

void SetupMover::restore_amide_hydrogen(pose::Pose& pose) {

	ReindexingMover forward_reindexer = ReindexingMover(1);
	ReindexingMover reverse_reindexer = ReindexingMover(-1);

	// Rotate the cycle forward.
	forward_reindexer.apply(pose);

	// Rebuild the amide hydrogen.
	PointPosition position;
	Size atom_index, residue_index = 2;
	conformation::Residue residue = pose.residue(residue_index);
	conformation::Conformation conformation = pose.conformation();
	id::NamedAtomID atom_name;

	atom_index = residue.atom_index("H");
	atom_name = id::NamedAtomID("H", residue_index);
	position = residue.build_atom_ideal(atom_index, conformation);
	pose.set_xyz(atom_name, position);

	// Rotate the cycle backward.
	reverse_reindexer.apply(pose);
}

} // End 'cycles' namespace.
} // End 'devel' namespace.


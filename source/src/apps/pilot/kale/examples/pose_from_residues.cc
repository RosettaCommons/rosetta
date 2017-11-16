// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>

static basic::Tracer TR( "apps.pilot.kale.examples" );

using namespace core;

int main(int argc, char **argv) {

	devel::init(argc, argv);

	pose::Pose source_pose, target_pose;
	Size source_index, target_index;
	Size residues, offset = 3;

	import_pose::pose_from_file(source_pose, "structures/marked_loop.8.pdb", core::import_pose::PDB_file);

	pose::remove_variant_type_from_pose_residue(
		source_pose, "LOWER_TERMINUS", 1);
	pose::remove_variant_type_from_pose_residue(
		source_pose, "UPPER_TERMINUS", source_pose.size());

	residues = source_pose.size();

	for ( Size i = 0; i < residues; i++ ) {
		source_index = 1 + (i - offset) % residues;
		target_index = 1 + i;

		conformation::Residue residue = source_pose.residue(source_index);
		target_pose.append_residue_by_bond(residue);
	}

	source_pose.dump_pdb("source.pdb");
	target_pose.dump_pdb("target.pdb");

}

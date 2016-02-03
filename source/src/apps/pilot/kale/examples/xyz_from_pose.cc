// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//! @example src/apps/pilot/kale/examples/xyz_from_pose.cc
//!
//! This example demonstrates how to determine the position of a specific atom 
//! within a pose.  The primary challenge is specifying the desired atom.

#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>

#include <numeric/xyzVector.hh>

// You need to include two header files in order to use the PointPosition 
// typedef:
//
//     #include <core/pose/Pose.hh>
//     #include <numeric/xyzVector.hh>
//
// The typedef itself is defined in <core/types.hh>, but this file is included 
// by <core/pose/Pose.hh>.  The underlying vector type behind this typedef is 
// defined in <numeric/xyzVector.hh>, so this file also needs to be included.  
// Without it, the compiler will generate an "incomplete aggregate type" error.

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kale.examples.XyzFromPose" );

int main(int argc, char* argv[]) {
		devel::init(argc, argv);

		pose::Pose pose;
		id::NamedAtomID atom;
		core::PointPosition position;

		import_pose::pose_from_file(pose, "structures/ubiquitin.pdb", core::import_pose::PDB_file);

		TR << pose.total_residue() << " residues loaded." << std::endl;

		atom = id::NamedAtomID("CA", 1);
		position = pose.xyz(atom);

		TR << "CA Position: (";
		TR << position.x() << ", ";
		TR << position.y() << ", ";
		TR << position.z() << ")" << std::endl;
}


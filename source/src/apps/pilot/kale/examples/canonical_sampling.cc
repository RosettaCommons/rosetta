// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

using namespace std;
using namespace protocols::canonical_sampling;

int main(int argc, char** argv) {
	devel::init(argc, argv);

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb(pose, "structures/linear/6.symmetry.pdb");

	//MetropolisHastingsMoverOP mover = new MetropolisHastingsMover;
}


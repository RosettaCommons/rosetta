// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <devel/cycles/SetupMover.hh>
#include <devel/cycles/ReindexingMover.hh>

#include <devel/init.hh>

using namespace core;
using namespace devel;

int main(int argc, char** argv) {

	devel::init(argc, argv);

	pose::Pose pose;
	cycles::SetupMover initializer;
	cycles::ReindexingMover reindexer(3);

	// If given `11.unmarked' instead of `8.marked', this program will output a 
	// structure with three missing residues.  More specifically, residues 1-3 
	// have the exact same coordinates as residues 6-8.  I'm not sure if this 
	// happens because the structure is unmarked or what.  This bug causes the 
	// Monte Carlo algorithm to explode after ~600 iterations.
	//
	// Before this code gets commited, I need to write rigorous tests for my 
	// cyclic movers.
	//
	// The pseudocode below might work a little better.  It basically hopes that 
	// the copy constructor is less fucked up than the assignment operator.  It's 
	// probably close to 2x slower, but that doesn't really matter.
	//
	// void apply(Pose &master_pose):
	// 		Pose source_pose (master_pose); // ie. copy ctor
	// 		Pose target_pose;
	//
	// 		// Build target_pose...
	//
	// 		master_pose.clear()
	// 		master_pose = target_pose;
	//
	// Here is some more information I've found regarding this bug:
	//
	// 1. The pseudocode above does not work any better.  Whatever the problem 
	// is, it aparently doesn't have to do with the assignment operator failing 
	// to properly copy the data structure.  Either that, or the copy constructor 
	// is fucking up in the same way.
	//
	// 2. The bug isn't due to the structure being unmarked.  The program fails 
	// in the same way when given `11.marked.pdb' as input.  Presumably, then, 
	// the bug is somehow brought on by the 3 extra residues.  
	//
	// 3. The `14.unmarked.pdb' structure also suffers this bug.  This isn't 
	// especially surprising, but it does indicate that the bug probably isn't 
	// due to some fluky feature of the `11.unmarked.pdb' structure.

	std::string input_pdb = "structures/cyclic/14.marked.pdb";
	import_pose::pose_from_file(pose, input_pdb, core::import_pose::PDB_file);

	initializer.apply(pose);
	reindexer.apply(pose);

	pose.dump_pdb("rotated_loop.14.pdb");
}


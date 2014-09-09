// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   debug_labontes_current_work.cc
/// @brief  This is simply a generic pilot app for testing things.
/// @author Labonte



// Package headers
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
//#include <core/pose/annotated_sequence.hh>
//#include <core/pose/PDBInfo.hh>
//#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/kinematics/FoldTree.hh>
//#include <core/kinematics/MoveMap.hh>
//#include <core/id/TorsionID.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>

//#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/simple_moves/PackRotamersMover.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <iostream>
//#include <algorithm>

// Construct random-number generator.
static numeric::random::RandomGenerator RG( 21 );  // the 6th triangular number


int main(int argc, char *argv[])
{
    try {
		using namespace std;
		using namespace core;
		using namespace scoring;
		using namespace import_pose;
		using namespace pose;

		// initialize core
		devel::init(argc, argv);

		// declare variables
		Pose pose, ref;

		// Make a test pose.
		//make_pose_from_sequence( pose, "AAAAAAAAAA", "fa_standard" );
		pose_from_pdb( ref, "/home/labonte/Workspace/Carbohydrates/MBP-G4_ref.pdb" );
		pose = ref;

		cout << pose << endl << endl;

		cout << "CA RMSD: " << CA_rmsd( pose, ref ) << endl;
		cout << "Heavy-atom ligand RMSD: " << non_peptide_heavy_atom_RMSD( pose, ref ) << endl;

		pose.set_phi( 372, 180.0 );

		cout << "CA RMSD: " << CA_rmsd( pose, ref ) << endl;
		cout << "Heavy-atom ligand RMSD: " << non_peptide_heavy_atom_RMSD( pose, ref ) << endl;

		//pose.dump_pdb( "" );
   } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}

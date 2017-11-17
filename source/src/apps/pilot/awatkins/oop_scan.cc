// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   oop_conformation.cc
/// @brief  Creates an OOP dimer and toys around with its dihedrals.
/// @author Watkins


// Package headers
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/io/carbohydrates/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
//#include <core/pose/PDBInfo.hh>
//#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/id/TorsionID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>

#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <iostream>
//#include <algorithm>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/ncbb/oop/OopCreatorMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

/*
code to scan OOP conformations--unfortunately, convergence is so bad that it isn't useful
Pose pose;

// Make an oop.
std::string sequence = "A[ALA:oop_pre]A[ALA:CtermProteinFull:oop_post]";
core::pose::make_pose_from_sequence( pose, sequence, core::chemical::FA_STANDARD );
std::cout << pose;

ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "mm_std" );

if ( scorefxn->has_zero_weight( core::scoring::atom_pair_constraint ) )
scorefxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );


std::cout << "Energies for pre angles with an optimized post" << std::endl << "\t";
for ( Real psi = -175; psi <= 180; psi += 10 ) {
std::cout << psi << "\t";
}
std::cout << std::endl;

for ( Real phi = -175; phi <= 180; phi += 10 ) {
std::cout << phi << "\t";
for ( Real psi = -175; psi <= 180; psi += 10 ) {

pose.set_phi( 1, phi );
pose.set_phi( 1, psi );

kinematics::MoveMapOP mm(new kinematics::MoveMap);
mm->set_bb( 2, true );

protocols::simple_moves::MinMover min( mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true );
min.apply( pose );
std::cout << (*scorefxn)(pose) << "\t";
}
std::cout << std::endl;
}
std::cout << std::endl;
*/

bool margin (
	core::Real val,
	core::Real comp,
	core::Real range
) {
	return ( ( val > comp - range && val < comp + range ) );
}

int
main( int argc, char *argv[] )
{
	try {
		using namespace std;
		using namespace utility;
		using namespace core;
		using namespace core::chemical;
		using namespace kinematics;
		using namespace scoring;
		using namespace import_pose;
		using namespace pose;
		using namespace protocols;
		using namespace simple_moves;
		using namespace oop;
		using namespace basic;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// initialize core
		devel::init( argc, argv );


		scoring::ScoreFunctionOP score_fxn = scoring::get_score_function();

		score_fxn->set_weight( core::scoring::atom_pair_constraint, 1 );
		score_fxn->set_weight( core::scoring::angle_constraint, 1.0 );
		score_fxn->set_weight( core::scoring::dihedral_constraint, 1 );//10.0 );

		Pose pose;
		std::string filename = option[in::file::s].value()[1];

		import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);

		bool anything_okay_for_orientation = false;
		Real mrgn = 30;

		while ( !anything_okay_for_orientation ) {

			for ( Size ii = 1; ii <= pose.size()-1; ++ii ) {
				if ( pose.residue( ii   ).type().name3() == "PRO" ) continue;
				if ( pose.residue( ii+1 ).type().name3() == "PRO" ) continue;
				if ( pose.residue( ii   ).type().has_variant_type( OOP_POST ) ) continue;

				// extrapolate where the OOP carbon positions would be
				if ( pose.residue( ii ).has( "H" ) && pose.residue( ii+1 ).has( "H" ) ) {
					Vector v1 = pose.residue( ii   ).atom( "H" ).xyz() - pose.residue( ii   ).atom( "N" ).xyz();
					Vector v2 = pose.residue( ii+1 ).atom( "H" ).xyz() - pose.residue( ii+1 ).atom( "N" ).xyz();
					if ( ( v2 - v1 ).length() < 1.6 ) {
						std::cout << ii << " and " << ( ii + 1 ) << " okay for distance" << std::endl;
						if ( margin( pose.phi( ii ), -145, mrgn ) && margin( pose.psi( ii ), -10, mrgn )
								&& ( margin( pose.phi( ii+1 ), -135, mrgn ) || margin( pose.phi( ii+1 ), 70, mrgn ) )
								&& margin( pose.psi( ii+1 ), 75, mrgn ) ) {
							std::cout << ii << " and " << ( ii + 1 ) << " okay for orientation" << std::endl;
							anything_okay_for_orientation = true;
						}
						if ( margin( pose.phi( ii ), 145, mrgn ) && margin( pose.psi( ii ), 10, mrgn )
								&& ( margin( pose.phi( ii+1 ), 135, mrgn ) || margin( pose.phi( ii+1 ), -70, mrgn ) )
								&& margin( pose.psi( ii+1 ), -75, mrgn ) ) {
							std::cout << ii << " and " << ( ii + 1 ) << " okay for orientation as Ds" << std::endl;
							anything_okay_for_orientation = true;
						}
					}
				}
			}
			mrgn += 1;
		}
		std::cout << "Final required margin was  " << mrgn << std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

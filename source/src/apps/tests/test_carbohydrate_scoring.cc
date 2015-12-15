// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    apps/tests/test_carbohydrate_scoring.cc
/// @brief   Application source code for testing carbohydrate scoring methods.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/MinMover.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>


using namespace std;
using namespace utility;
using namespace core;
using namespace chemical::carbohydrates;
using namespace pose;
using namespace import_pose;
using namespace scoring;
using namespace protocols;


string const PATH = "input/";


Real
output_score( Pose & sugar, core::uint res_num, ScoreFunction const & sf )
{
	cout << " Phi/Psi: " << sugar.phi( res_num ) << '/' << sugar.psi( res_num );
	cout << "  Total Score: " << sf( sugar );
	cout << "  Sugar BB Score (res " << res_num << "): ";
	Real const sugar_bb_score( sugar.energies().residue_total_energies( res_num )[ sugar_bb ] );
	cout << sugar_bb_score << endl;

	return sugar_bb_score;
}

void
sample_torsions( Pose & pose, core::uint res_num, ScoreFunction const & sf )
{
	CarbohydrateInfoCOP info( pose.residue( res_num ).carbohydrate_info() );

	cout << "Residue " << res_num << " is ";
	if ( info->is_alpha_sugar() ) {
		cout << "an alpha sugar." << endl;
	} else {
		cout << "a beta sugar." << endl;
	}
	cout << "It is attached to:" << endl;
	core::uint const parent( carbohydrates::find_seqpos_of_saccharides_parent_residue( pose.residue( res_num ) ) );
	cout << "Residue " << parent << ": " << pose.residue( parent ).name() << endl;
	cout << "Sampling glycosidic bonds for residue " << res_num << "..." << endl;

	Real worst_sugar_bb_score( 0.0 ), best_sugar_bb_score( 9999.0 );  // arbitrarily large number
	Angle worst_phi, worst_psi, best_phi, best_psi;

	for ( Angle phi( -180.0 ); phi <= 180.0; phi += 30.0 ) {
		pose.set_phi( res_num, phi );

		for ( Angle psi( -180.0 ); psi <= 180.0; psi += 30.0 ) {
			pose.set_psi( res_num, psi );

			Real const sugar_bb_score( output_score( pose, res_num, sf ) );
			if ( sugar_bb_score < best_sugar_bb_score ) {
				best_sugar_bb_score = sugar_bb_score;
				best_phi = phi;
				best_psi = psi;
			}
			if ( sugar_bb_score > worst_sugar_bb_score ) {
				worst_sugar_bb_score = sugar_bb_score;
				worst_phi = phi;
				worst_psi = psi;
			}
		}
		cout << endl;
	}
	cout << "Best Phi/Psi combination: " << best_phi << '/' << best_psi;
	cout << " Sugar BB Score: " << best_sugar_bb_score << endl;
	cout << "Worst Phi/Psi combination: " << worst_phi << '/' << worst_psi;
	cout << " Sugar BB Score: " << worst_sugar_bb_score << endl;
}


int
main( int argc, char *argv[] )
{
	try {
		// Initialize core.
		devel::init( argc, argv );

		// Declare Pose variables.
		Pose maltotriose, lactose, LeX, O_linked_glycan;

		// Set up ScoreFunctions.
		ScoreFunctionOP sf_full( get_score_function() );
		ScoreFunctionOP sf_bb( new ScoreFunction );
		sf_bb->set_weight( sugar_bb, 1.0 );

		// Set up MinMovers.
		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_nu( false );
		mm->set_branches( true );
		simple_moves::MinMover minimizer( mm, sf_full, "dfpmin", 0.01, true );
		cout << minimizer << endl;


		// Begin testing. /////////////////////////////////////////////////////

		cout << "Importing maltotriose (D-alpha1->4eq linkage):" << endl;

		pose_from_pdb( maltotriose, PATH + "maltotriose.pdb" );

		sample_torsions( maltotriose, 2, *sf_full );
		cout << "(The best angles should be close to 60/90.)" << endl;
		cout << "(The worst angles should be close to -90/-120 or 0.)" << endl;

		cout << "Setting Phi/Psi of residue 2 to the best angles...." << endl;
		maltotriose.set_phi( 2, 60.0 );
		maltotriose.set_psi( 2, 90.0 );

		sample_torsions( maltotriose, 3, *sf_full );
		cout << "(The best angles should be close to 60/90.)" << endl;
		cout << "(The worst angles should be close to -90/-120 or 0.)" << endl;


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing lactose (D-beta1->4eq linkage):" << endl;

		pose_from_pdb( lactose, PATH + "lactose.pdb" );

		sample_torsions( lactose, 2, *sf_full );
		cout << "(The best angles should be close to -60/90.)" << endl;
		cout << "(The worst angles should be close to 120/-120 or 0.)" << endl;

		cout << "Setting lactose's glycosidic bond on slope toward minimum..." << endl;
		lactose.set_phi( 2, -120.0 );
		lactose.set_psi( 2, 30.0 );

		output_score( lactose, 2, *sf_full );

		cout << "Minimizing lactose..." << endl;

		minimizer.apply( lactose );

		output_score( lactose, 2, *sf_full );

		cout << "Repeating process using the sugar_bb scoring term only..." << endl;

		minimizer.score_function( sf_bb );

		lactose.set_phi( 2, -120.0 );
		lactose.set_psi( 2, 30.0 );

		output_score( lactose, 2, *sf_bb );

		minimizer.apply( lactose );

		output_score( lactose, 2, *sf_bb );

		cout << "(The minimized angles should be close to -60/90.)" << endl;


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing LewisX (D-beta1->4eq main-chain linkage and L-alpha1->3eq branch linkage):" << endl;

		pose_from_pdb( LeX, PATH + "Lex.pdb" );

		output_score( LeX, 3, *sf_full );

		cout << "Setting LewisX's glycosidic bond on slope toward minimum..." << endl;
		LeX.set_phi( 3, -120.0 );
		LeX.set_psi( 3, -30.0 );

		output_score( LeX, 3, *sf_full );

		cout << "Minimizing LewisX..." << endl;

		minimizer.score_function( sf_full );
		minimizer.apply( LeX );

		output_score( LeX, 3, *sf_full );

		cout << "Repeating process using the sugar_bb scoring term only..." << endl;

		minimizer.score_function( sf_bb );

		LeX.set_phi( 3, -120.0 );
		LeX.set_psi( 3, -30.0 );

		output_score( LeX, 3, *sf_bb );

		minimizer.apply( LeX );

		output_score( LeX, 3, *sf_bb );

		cout << "(The minimized angles should be close to -60/-90.)" << endl;


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing sample O-linked glycan (D-alpha1->OSer linkage):" << endl;

		pose_from_pdb( O_linked_glycan, PATH + "O_glycan.pdb" );

		output_score( O_linked_glycan, 4, *sf_full );

		cout << "Setting sample O-linked glycan's glycosidic bond on slope toward minimum..." << endl;
		O_linked_glycan.set_phi( 4, -120.0 );

		output_score( O_linked_glycan, 4, *sf_full );

		cout << "Minimizing O-linked glycan..." << endl;

		minimizer.score_function( sf_full );
		minimizer.apply( O_linked_glycan );

		output_score( O_linked_glycan, 4, *sf_full );

		cout << "Repeating process using the sugar_bb scoring term only..." << endl;

		minimizer.score_function( sf_bb );

		O_linked_glycan.set_phi( 4, -120.0 );

		output_score( O_linked_glycan, 4, *sf_bb );

		minimizer.apply( O_linked_glycan );

		output_score( O_linked_glycan, 4, *sf_bb );

		cout << "(The minimized angles should be close to 60/NA.)" << endl;

	} catch ( excn::EXCN_Base const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    apps/tests/test_carbohydrate_scoring.cc
/// @brief   Application source code for testing carbohydrate scoring methods.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
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
using namespace core::import_pose;
using namespace chemical::carbohydrates;
using namespace pose;
using namespace pose::carbohydrates;
using namespace import_pose;
using namespace scoring;
using namespace protocols;


string const PATH = "input/";


Real
output_score( Pose & sugar, core::uint res_num, ScoreFunction const & sf, bool silently = false )
{
	Real const total_score( sf( sugar ) );
	Real const sugar_bb_score( sugar.energies().residue_total_energies( res_num )[ sugar_bb ] );

	if ( ! silently ) {
		if (  sugar.glycan_tree_set()->has_exocyclic_glycosidic_linkage( res_num ) ) {
			cout << " Phi/Psi/Omega: " <<
				sugar.phi( res_num ) << '/' << sugar.psi( res_num ) << '/' << sugar.omega( res_num );
		} else {
			cout << " Phi/Psi: " << sugar.phi( res_num ) << '/' << sugar.psi( res_num );
		}
		cout << "  Total Score: " << total_score;
		cout << "  Sugar BB Score (res " << res_num << "): " << sugar_bb_score << endl;
	}

	return total_score;
}

void
sample_torsions( Pose & pose, core::uint res_num, ScoreFunction const & sf )
{
	bool const silently( true );
	CarbohydrateInfoCOP info( pose.residue( res_num ).carbohydrate_info() );

	cout << "Residue " << res_num << " is ";
	if ( info->is_alpha_sugar() ) {
		cout << "an alpha sugar." << endl;
	} else {
		cout << "a beta sugar." << endl;
	}
	cout << "It is attached to:" << endl;
	core::uint const parent( pose.glycan_tree_set()->get_parent( res_num ) ) ;
	cout << "Residue " << parent << ": " << pose.residue( parent ).name() << endl;
	cout << "Sampling glycosidic bonds for residue " << res_num << "..." << endl;

	Real worst_score( 0.0 ), best_score( 9999.0 );  // arbitrarily large number
	Angle worst_phi(-999), worst_psi(-999), worst_omega(-999), best_phi(-999), best_psi(-999), best_omega(-999);
	bool const sample_omega( pose.glycan_tree_set()->has_exocyclic_glycosidic_linkage( res_num ) );

	for ( Angle phi( -180.0 ); phi <= 180.0; phi += 15.0 ) {
		pose.set_phi( res_num, phi );

		for ( Angle psi( -180.0 ); psi <= 180.0; psi += 15.0 ) {
			pose.set_psi( res_num, psi );

			if ( sample_omega ) {
				for ( Angle omega( -180.0 ); omega <= 180.0; omega += 15.0 ) {
					pose.set_omega( res_num, omega );

					Real const total_score( output_score( pose, res_num, sf, silently ) );
					if ( total_score < best_score ) {
						best_score = total_score;
						best_phi = phi;
						best_psi = psi;
						best_omega = omega;
					}
					if ( total_score > worst_score ) {
						worst_score = total_score;
						worst_phi = phi;
						worst_psi = psi;
						worst_omega = omega;
					}
				}
			} else {
				Real const total_score( output_score( pose, res_num, sf, silently ) );
				if ( total_score < best_score ) {
					best_score = total_score;
					best_phi = phi;
					best_psi = psi;
				}
				if ( total_score > worst_score ) {
					worst_score = total_score;
					worst_phi = phi;
					worst_psi = psi;
				}
			}
		}
		//cout << endl;
	}
	if ( sample_omega ) {
		cout << "Best Phi/Psi/Omega combination: " << best_phi << '/' << best_psi << '/' << best_omega;
	} else {
		cout << "Best Phi/Psi combination: " << best_phi << '/' << best_psi;
	}
	cout << " Total Score: " << best_score << endl;
	if ( sample_omega ) {
		cout << "Worst Phi/Psi/Omega combination: " << worst_phi << '/' << worst_psi << '/' << worst_omega;
	} else {
		cout << "Worst Phi/Psi combination: " << worst_phi << '/' << worst_psi;
	}
	cout << " Total Score: " << worst_score << endl;
}


int
main( int argc, char *argv[] )
{
	try {
		// Initialize core.
		devel::init( argc, argv );

		// Declare Pose variables.
		Pose maltobiose, lactose, isomaltose, LeX, bisected_man, N_linked_glycan, O_linked_glycan;

		// Set up ScoreFunctions.
		ScoreFunctionOP sf_full( get_score_function() );
		ScoreFunctionOP sf_no_bb( sf_full->clone() );
		sf_no_bb->set_weight( sugar_bb, 0.0 );
		ScoreFunctionOP sf_bb( new ScoreFunction );
		sf_bb->set_weight( sugar_bb, 1.0 );

		// Set up MinMovers.
		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_nu( false );
		mm->set_branches( true );
		simple_moves::MinMover minimizer;
		minimizer.movemap( mm );
		minimizer.score_function( sf_bb );
		//minimizer.min_type( "linmin_iterated" );
		minimizer.tolerance( 0.001 );


		// Begin testing. /////////////////////////////////////////////////////

		// Test Scoring of full range of angles.

		cout << "Importing maltotriose (D-alpha1->4eq linkage):" << endl << endl;

		pose_from_file( maltobiose, PATH + "maltobiose.pdb", PDB_file );

		cout << "Sampling WITHOUT the sugar_bb scoring term..." << endl;

		sample_torsions( maltobiose, 2, *sf_no_bb );
		cout << "(The best angles should be close to 75/90.)" << endl;
		cout << "(The worst angles should be close to -105/-120 or -15.)" << endl;

		cout << endl;
		cout << "Repeating sampling WITH the sugar_bb term..." << endl;

		sample_torsions( maltobiose, 2, *sf_full );
		cout << "(The best angles should be close to 75/90.)" << endl;
		cout << "(The worst angles should be close to -105/-120 or -15.)" << endl;

		cout << endl;
		cout << "Repeating sampling with the sugar_bb term ONLY..." << endl;

		sample_torsions( maltobiose, 2, *sf_bb );
		cout << "(The best angles should be close to 75/90.)" << endl;
		cout << "(The worst angles should be close to -105/-120 or -15.)" << endl;


		cout << endl << "-------------------------------------------------------------------------------" << endl;

		cout << "Importing lactose (D-beta1->4eq linkage):" << endl << endl;

		pose_from_file( lactose, PATH + "lactose.pdb", PDB_file );

		sample_torsions( lactose, 2, *sf_bb );
		cout << "(The best angles should be close to -60/90.)" << endl;
		cout << "(The worst angles should be close to 120/-120 or -15.)" << endl;

		cout << endl << "-------------------------------------------------------------------------------" << endl;

		cout << "Importing isomaltose (D-alpha1->6 linkage):" << endl << endl;

		pose_from_file( isomaltose, PATH + "isomaltose.pdb", PDB_file );

		sample_torsions( isomaltose, 2, *sf_bb );
		cout << "(The best angles should be close to 75/180/-60.)" << endl;
		cout << "(The worst angles should be close to -105/0/-120.)" << endl;

		cout << endl << "-------------------------------------------------------------------------------" << endl;


		// Test minimization.

		cout << "Setting lactose's glycosidic bond on slope toward minimum..." << endl;
		lactose.set_phi( 2, -120.0 );
		lactose.set_psi( 2, 30.0 );

		output_score( lactose, 2, *sf_full );

		cout << "Minimizing lactose using sugar_bb scoring term only..." << endl;

		minimizer.apply( lactose );

		output_score( lactose, 2, *sf_full );

		cout << "(The minimized angles should be close to -60/90.)" << endl;

		cout << endl << "-------------------------------------------------------------------------------" << endl;

		cout << "Setting isomaltose's glycosidic bond on slope toward minimum..." << endl;
		isomaltose.set_phi( 2, 30.0 );
		isomaltose.set_psi( 2, -150.0 );
		isomaltose.set_omega( 2, -30.0 );

		output_score( isomaltose, 2, *sf_full );

		cout << "Minimizing isomaltose using sugar_bb scoring term only..." << endl;

		minimizer.apply( isomaltose );

		output_score( isomaltose, 2, *sf_full );

		cout << "(The minimized angles should be close to 75/180/-60.)" << endl;

		cout << endl << "-------------------------------------------------------------------------------" << endl;

		cout << "Importing LewisX (D-beta1->4eq main-chain linkage and L-alpha1->3eq branch linkage):" << endl << endl;

		pose_from_file( LeX, PATH + "Lex.pdb", PDB_file );

		cout << "Setting LewisX's glycosidic bonds on slope toward minimum..." << endl;
		LeX.set_phi( 2, -120.0 );
		LeX.set_psi( 2, 30.0 );
		LeX.set_phi( 3, -120.0 );
		LeX.set_psi( 3, -30.0 );

		output_score( LeX, 2, *sf_full );
		output_score( LeX, 3, *sf_full );

		cout << "Minimizing LewisX using sugar_bb scoring term only..." << endl;

		minimizer.apply( LeX );

		output_score( LeX, 2, *sf_full );
		output_score( LeX, 3, *sf_full );

		cout << "(The minimized angles should be close to -60/90 for residue 2.)" << endl;
		cout << "(The minimized angles should be close to -75/-90 for residue 3.)" << endl;

		cout << endl << "-------------------------------------------------------------------------------" << endl;

		cout << "Creating oligosaccharide with multiple branches off a single residue, ";
		cout << "(D-alpha1->3eq main-chain linkage and D-beta1->4eq and D-alpha1->6 branch linkages):" << endl << endl;

		make_pose_from_saccharide_sequence( bisected_man,
			"a-D-Manp-(1->3)-[a-D-Manp-(1->6)]-[b-d-GlcpNAc-(1->4)]-b-D-Manp" );

		cout << "Setting bisected oligosaccharide's bonds on slope toward minimum..." << endl;
		bisected_man.set_phi( 2, 30.0 );
		bisected_man.set_psi( 2, -30.0 );
		bisected_man.set_phi( 3, -30.0 );
		bisected_man.set_psi( 3, 30.0 );
		bisected_man.set_phi( 4, 30.0 );
		bisected_man.set_psi( 4, -150.0 );
		bisected_man.set_omega( 4, -90.0 );

		output_score( bisected_man, 2, *sf_full );
		output_score( bisected_man, 3, *sf_full );
		output_score( bisected_man, 4, *sf_full );

		cout << "Minimizing using sugar_bb scoring term only..." << endl;

		minimizer.apply( bisected_man );

		output_score( bisected_man, 2, *sf_full );
		output_score( bisected_man, 3, *sf_full );
		output_score( bisected_man, 4, *sf_full );

		cout << "(The minimized angles should be close to 75/-90 for residue 2.)" << endl;
		cout << "(The minimized angles should be close to -60/90 for residue 3.)" << endl;
		cout << "(The minimized angles should be close to 75/180/-60 for residue 4.)" << endl;

		cout << endl << "-------------------------------------------------------------------------------" << endl;

		cout << "Importing sample N-linked glycan (D-beta1->NAsn linkage):" << endl << endl;

		pose_from_file( N_linked_glycan, PATH + "N-linked_14-mer_glycan.pdb", PDB_file );

		cout << "Setting sample N-linked glycan's glycosidic bond on slope toward minimum..." << endl;
		N_linked_glycan.set_phi( 6, -90.0 );

		output_score( N_linked_glycan, 6, *sf_full );

		cout << "Minimizing O-linked glycan using sugar_bb scoring term only..." << endl;

		minimizer.apply( N_linked_glycan );

		output_score( N_linked_glycan, 6, *sf_full );

		cout << "(The minimized angles should be close to -60/NA.)" << endl;

		cout << endl << "-------------------------------------------------------------------------------" << endl;

		cout << "Importing sample O-linked glycan (D-alpha1->OSer linkage):" << endl << endl;

		pose_from_file( O_linked_glycan, PATH + "O_glycan.pdb", PDB_file );

		cout << "Setting sample O-linked glycan's glycosidic bond on slope toward minimum..." << endl;
		O_linked_glycan.set_phi( 4, -120.0 );

		output_score( O_linked_glycan, 4, *sf_full );

		cout << "Minimizing O-linked glycan using sugar_bb scoring term only..." << endl;

		minimizer.apply( O_linked_glycan );

		output_score( O_linked_glycan, 4, *sf_full );

		cout << "(The minimized angles should be close to 75/NA.)" << endl;

	} catch (excn::Exception const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
// Protocol Headers
#include <devel/init.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <protocols/pockets/FingerprintMultifunc.hh>

#include <basic/options/option_macros.hh>
#include <utility/excn/Exceptions.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>

#include <basic/options/option_macros.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;


OPT_KEY( String, input_ligand_file1 )
OPT_KEY( String, input_ligand_file2 )
OPT_KEY( Real, phi_increment )
OPT_KEY( Real, psi_increment )

int main( int argc, char * argv [] ) {
	try{

	NEW_OPT( input_ligand_file1, "ligand1 file name", "ligand1.pdb" );
	NEW_OPT( input_ligand_file2, "ligand2 file name", "ligand2.pdb" );
	NEW_OPT( phi_increment, "phi increment", 2.0 );
	NEW_OPT( psi_increment, "psi increment", 2.0 );

  devel::init(argc, argv);

	std::string const input_ligand1 = option[ input_ligand_file1 ];
	std::string const input_ligand2 = option[ input_ligand_file2 ];
	core::Real const phi_inc  = option[ phi_increment ];
	core::Real const psi_inc  = option[ psi_increment ];

	/*
	protocols::pockets::NonPlaidFingerprint npf;
	npf.setup_from_PlaidFingerprint(pf1);
	std::string np_output_pdbname = "npf_output_pdb.pdb";
	std::string np_output_filename = "npf_output_file.txt";
	npf.print_to_pdb(np_output_pdbname);
	npf.print_to_file(np_output_filename);

	pose::Pose small_mol_pose1;
<<<<<<< HEAD
	core::import_pose::pose_from_pdb( small_mol_pose1, input_ligand1 );
	core::grid::Pockets::PlaidFingerprint pf1( small_mol_pose1, npf );
=======
	io::pdb::pose_from_pdb( small_mol_pose1, input_ligand1 );
	protocols::pockets::PlaidFingerprint pf1( small_mol_pose1, npf );
>>>>>>> origin/master
	//	pf1.build_from_pose( small_mol_pose1 );
	std::string p_output_pdbname = "pf1_output_pdb.pdb";
	pf1.print_to_pdb(p_output_pdbname);

	pose::Pose small_mol_pose2;
<<<<<<< HEAD
	core::import_pose::pose_from_pdb( small_mol_pose2, input_ligand2 );
	core::grid::Pockets::PlaidFingerprint pf2( small_mol_pose2, npf );
=======
	io::pdb::pose_from_pdb( small_mol_pose2, input_ligand2 );

	protocols::pockets::PlaidFingerprint pf2( small_mol_pose2, npf );
>>>>>>> origin/master
	//pf2.build_from_pose( small_mol_pose2 );

	std::string	pf2_output_filename = "pf2_output_file.txt";
	std::string pf2_output_pdbname = "pf2_output_pdb.pdb";
	pf2.print_to_file(pf2_output_filename);
	pf2.print_to_pdb(pf2_output_pdbname);

	core::Real unaligned_score = pf2.fp_compare( npf );
	std::cout << "Unaligned score is: " << unaligned_score << std::endl;

 //Move origin
	numeric::xyzVector<core::Real> move_npf;
	move_npf.x() = 42.48;
	move_npf.y() = 11.60;
	move_npf.z() = 29.54;
	npf.move_origin(move_npf);
	np_output_filename = "after_move_file.txt";
	np_output_pdbname = "after_move_pdb.pdb";
	npf.print_to_file(np_output_filename);
	npf.print_to_pdb(np_output_pdbname);
	std::cout << "After move origin score is: " <<  pf.fp_compare( npf ) << std::endl;


	core::Real optimal_score, optimal_phi, optimal_psi;
	optimal_score =  pfi.find_optimal_rotation( npf, phi_inc, psi_inc, optimal_phi, optimal_psi );
	std::cout << "Optimal_score: " << optimal_score << " Optimal_angles: " << optimal_phi <<" " << optimal_psi <<std::endl;

	//Apply rotation
	npf.apply_rotation(optimal_phi,optimal_psi);
	std::cout << "apply_rotation score : " <<  pfi.fp_compare( npf ) << std::endl;
	np_output_filename = "rotated_file.txt";
	np_output_pdbname = "rotated_pdb.pdb";
	npf.print_to_file(np_output_filename);
	npf.print_to_pdb(np_output_pdbname);
	*/

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl; 
	} 
	return 0;

}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test_C-terminal_conjugation.cc
/// @brief  This is an integration test app for testing reverse-fold-tree generation code.
/// @author Labonte <JWLabonte@jhu.edu>

// includes
#include <devel/init.hh>

//#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>

#include <utility/excn/Exceptions.hh>

#include <iostream>


using namespace std;
using namespace core;
using namespace import_pose;
using namespace pose;


string const INPATH = "input/";
string const OUTPATH = "output/";


int
main( int argc, char *argv[] )
{
	using namespace std;

	try {
		// Initialize core.
		devel::init( argc, argv );

		// Import test poses.
		Pose pose, O_linked_glycan_pose, O_linked_centroid_glycan_pose, conjugated_pose1, conjugated_pose2, conjugated_pose3;
		pose_from_file( pose, INPATH + "test_normal1.pdb" , core::import_pose::PDB_file);
		pose_from_file( conjugated_pose1, INPATH + "test_c_conj1.pdb" , core::import_pose::PDB_file);
		pose_from_file( conjugated_pose2, INPATH + "UBQ_E2_0001.pdb" , core::import_pose::PDB_file);
		pose_from_file( conjugated_pose3, INPATH + "UBQ_E2_0003.pdb" , core::import_pose::PDB_file);
		pose_from_file( O_linked_glycan_pose, INPATH + "test_O_glycan.pdb" , core::import_pose::PDB_file);
		chemical::ResidueTypeSetCOP centroid_rts = chemical::ChemicalManager::get_instance()->residue_type_set(chemical::CENTROID);
		pose_from_file( O_linked_centroid_glycan_pose, *centroid_rts, INPATH + "test_O_glycan_centroid.pdb" , core::import_pose::PDB_file);

		cout << "PDB file without LINK record" << endl;
		cout << pose << endl;
		cout << "1: " << pose.residue( 1 ).name() << endl;
		cout << "84: " << pose.residue( 84 ).name() << endl;
		cout << "146: " << pose.residue( 146 ).name() << endl;
		cout << "147: " << pose.residue( 147 ).name() << endl;
		cout << "203: " << pose.residue( 203 ).name() << endl;
		cout << "204: " << pose.residue( 204 ).name() << endl;
		cout << "280: " << pose.residue( 280 ).name() << endl;

		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Same PDB file with LINK record (O-linked)" << endl;
		cout << "" << conjugated_pose1 << endl;
		cout << "1: " << conjugated_pose1.residue( 1 ).name() << endl;
		cout << "84: " << conjugated_pose1.residue( 84 ).name() << endl;
		cout << "146: " << conjugated_pose1.residue( 146 ).name() << endl;
		cout << "147: " << conjugated_pose1.residue( 147 ).name() << endl;
		cout << "203: " << conjugated_pose1.residue( 203 ).name() << endl;
		cout << "204: " << conjugated_pose1.residue( 204 ).name() << endl;
		cout << "280: " << conjugated_pose1.residue( 280 ).name() << endl;

		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "S-linked example" << endl;
		cout << "" << conjugated_pose2 << endl;

		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Another S-linked example" << endl;
		cout << "" << conjugated_pose3 << endl;

		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "An O-linked glycopeptide" << endl;
		cout << "" << O_linked_glycan_pose << endl;
		cout << "1: " << O_linked_glycan_pose.residue( 1 ).name() << endl;
		cout << "2: " << O_linked_glycan_pose.residue( 2 ).name() << endl;
		cout << "3: " << O_linked_glycan_pose.residue( 3 ).name() << endl;
		cout << "4: " << O_linked_glycan_pose.residue( 4 ).name() << endl;

		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "A centroid O-linked glycopeptide" << endl;
		cout << "" << O_linked_centroid_glycan_pose << endl;
		cout << "1: " << O_linked_centroid_glycan_pose.residue( 1 ).name() << endl;
		cout << "2: " << O_linked_centroid_glycan_pose.residue( 2 ).name() << endl;
		cout << "3: " << O_linked_centroid_glycan_pose.residue( 3 ).name() << endl;
		cout << "4: " << O_linked_centroid_glycan_pose.residue( 4 ).name() << endl;

		conjugated_pose1.dump_pdb( OUTPATH + "O-linked_Ubi_peptide.pdb" );
		conjugated_pose3.dump_pdb( OUTPATH + "S-linked_Ubi_peptide.pdb" );
		O_linked_glycan_pose.dump_pdb( OUTPATH + "O-linked_glycopeptide.pdb" );
		O_linked_centroid_glycan_pose.dump_pdb( OUTPATH + "O-linked_centroid_glycopeptide.pdb" );
	} catch (utility::excn::Exception const & e ) {
		cerr << "caught exception " << e.msg() << endl;
		return -1;
	}
	return 0;
}

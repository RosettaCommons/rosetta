// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test_CarbohydrateInfo.cc
/// @brief   Pilot application source code for testing CarbohydrateInfo.
/// @author  labonte


// Project headers
#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/RingConformerSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>

using namespace std;
using namespace core;
using namespace pose;
using namespace import_pose;
using namespace chemical;
using namespace conformation;


string const PATH = "input/";


void
test_sugar(Pose & sugar)
{
	cout << endl << sugar << endl;

	cout << "Sequences:" << endl;
	Size n_chains = sugar.conformation().num_chains();
	for (core::uint i = 1; i <= n_chains; ++i) {
		cout << " Chain " << i << ": ";
		cout << sugar.chain_sequence(i) << endl;
	}

	cout << endl << "Residue Info:" << endl;
	Size n_res = sugar.total_residue();
	for (core::uint i = 1; i <= n_res; ++i) {
		Residue res = sugar.residue(i);
		cout << "PDB ID: " << sugar.pdb_info()->pose2pdb(i) << ": " << res << endl << endl;
	}
}


int
main(int argc, char *argv[])
{
	try {
		// Initialize core.
		devel::init(argc, argv);

		// Declare variables.
		Pose maltotriose, isomaltose, lactose, amylopectin, glycopeptide, glucosamine, N_linked_14_mer, psicose,
				neuraminate, Lex, GalCer;

		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing maltotriose:" << endl;

		pose_from_pdb(maltotriose, PATH + "maltotriose.pdb");

		test_sugar(maltotriose);


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing isomaltose:" << endl;

		pose_from_pdb(isomaltose, PATH + "isomaltose.pdb");

		test_sugar(isomaltose);


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing lactose:" << endl;

		pose_from_pdb(lactose, PATH + "lactose.pdb");

		test_sugar(lactose);


		/*cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Creating maltotriose from sequence:" << endl;

		ResidueTypeSetCAP residue_set(ChemicalManager::get_instance()->residue_type_set("fa_standard"));
		make_pose_from_saccharide_sequence(maltotriose, "alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp", *residue_set);

		// TODO: Figure out what is wrong from sugar poses made from sequences.
		//test_sugar(maltotriose);
		cout << endl << maltotriose << endl;

		cout << "Sequences:" << endl;
		cout << " Chain 1: ";
		cout << maltotriose.chain_sequence(1) << endl;*/


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing branched amylopectin fragment:" << endl;

		pose_from_pdb(amylopectin, PATH + "amylopectin_fragment.pdb");

		test_sugar(amylopectin);


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing N-glycosylated sample:" << endl;

		pose_from_pdb(glycopeptide, PATH + "glycosylated_peptide.pdb");

		test_sugar(glycopeptide);
		cout << endl << glycopeptide << endl;

		cout << "Sequences:" << endl;
		for (core::uint i = 1; i <= 2; ++i) {
			cout << " Chain " << i << ": ";
			cout << glycopeptide.chain_sequence(i) << endl;
		}


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing glucosamine:" << endl;

		pose_from_pdb(glucosamine, PATH + "GlcN.pdb");

		test_sugar(glucosamine);

		cout << *glucosamine.residue(1).type().ring_conformer_set();


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing N-glycosylated 14-mer:" << endl;

		pose_from_pdb(N_linked_14_mer, PATH + "N-linked_14-mer_glycan.pdb");

		test_sugar(N_linked_14_mer);


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing beta-D-psicopyranose:" << endl;

		pose_from_pdb(psicose, PATH + "beta-psicose.pdb");

		test_sugar(psicose);

		cout << *psicose.residue(1).type().ring_conformer_set();


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing alpha-neuraminic acid:" << endl;

		pose_from_pdb(neuraminate, PATH + "Neu.pdb");

		test_sugar(neuraminate);


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing Lewisx:" << endl;

		pose_from_pdb(Lex, PATH + "Lex.pdb");

		test_sugar(Lex);


		cout << "---------------------------------------------------------------------------------------------" << endl;
		cout << "Importing GalCer:" << endl;

		pose_from_pdb(GalCer, PATH + "GalCer.pdb");

		test_sugar(GalCer);

	} catch (utility::excn::EXCN_Base const & e) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}

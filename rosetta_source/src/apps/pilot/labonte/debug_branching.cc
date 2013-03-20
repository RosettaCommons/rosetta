// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    debug_branching.cc
/// @brief   App for figuring out how to do branching.
/// @author  labonte


// Project headers
#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>
#include <iostream>


using namespace std;
using namespace utility;
using namespace core;
using namespace chemical;
using namespace conformation;
using namespace pose;
using namespace import_pose;


string const PATH = "/home/labonte/Workspace/";


int
main(int argc, char *argv[])
{
    try {
	// Initialize core.
	devel::init(argc, argv);

	vector1<string> params_paths;
	Pose pose;

	params_paths.push_back(PATH + "Carbohydrates/LG1.params");
	params_paths.push_back(PATH + "Carbohydrates/LG2.params");

	ResidueTypeSetOP type_set(ChemicalManager::get_instance()->nonconst_residue_type_set("fa_standard"));
	type_set->read_files(params_paths,
			ChemicalManager::get_instance()->atom_type_set("fa_standard"),
			ChemicalManager::get_instance()->element_set("fa_standard"),
			ChemicalManager::get_instance()->mm_atom_type_set("fa_standard"),
			ChemicalManager::get_instance()->orbital_type_set("fa_standard"));

	pose_from_pdb(pose, *type_set, PATH + "Carbohydrates/LG_0001.pdb");

	cout << pose << endl << endl;

	Size n_residues = pose.total_residue();
	for (core::uint res_num = 1; res_num <= n_residues; ++res_num) {
		ResidueCAP residue = & pose.residue(res_num);
		cout << *residue << endl << endl;
	}

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}

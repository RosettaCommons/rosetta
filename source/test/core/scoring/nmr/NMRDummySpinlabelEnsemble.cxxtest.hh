// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDummySpinlabelEnsemble.cxxtest.hh
/// @brief   unit test for class NMRDummySpinlabelEnsemble that stores and handles data of
///          a virtual spinlabel. This is represented as a collection of N conformers with its
///          heavy atom coordinates and number of observations stored in the database. The
///          ensemble has a set of parameters associated with it to perform a clash score
///          calculation which is based on the measurement of atom pair distances to the
///          neighboring amino acids.
/// @details Last modified: 12/04/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>
#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <numeric>

static basic::Tracer TR("core.scoring.nmr.NMRDummySpinlabelEnsemble.cxxtest");

class NMRDummySpinlabelEnsembleTests : public CxxTest::TestSuite {

private:
	core::pose::Pose cam_;
	std::string ensemble_conformers_database_file_;

public:
	/// @brief Setup Test
	void setUp() {

		// Initialize core & options system
		core_init();

		// Load pose from pdb
		core::import_pose::pose_from_file(cam_, "core/scoring/nmr/1cdl.pdb", core::import_pose::PDB_file);
		ensemble_conformers_database_file_ = basic::database::full_name("chemical/residue_type_sets/fa_standard/residue_types/spin_labels/R1A_rotamers_nonclashing_moe_saved.pdb");
	}

	void tearDown() {
		cam_.clear();
	}

	/// @brief test creation of NMRDummySpinlabelEnsemble and NMRSpinlabel objects
	///        and access of NMRDummySpinlabelEnsemble data
	void test_instantiation_and_access() {
		using namespace core::scoring::nmr;

		TR << "Testing instantiation of NMRDummySpinlabelEnsemble." << std::endl;

		NMRDummySpinlabelEnsembleOP r1a_virt_ensemble;
		core::chemical::ResidueTypeOP restype( new core::chemical::ResidueType( *((core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard"))->get_representative_type_name3("R1A")) ) );
		try {
			r1a_virt_ensemble = NMRDummySpinlabelEnsembleOP( new NMRDummySpinlabelEnsemble(ensemble_conformers_database_file_, *restype) );
		} catch ( utility::excn::Exception & excn ) {
			TR << "Caught an exception during creation of NMRDummySpinlabelEnsemble: " << excn.msg() << std::endl;
			utility_exit_with_message("Exiting program ...");
		} catch (...) {
			utility_exit_with_message("Caught an exception during creation of NMRDummySpinlabelEnsemble. Exiting program ...");
		}
		TS_ASSERT_EQUALS(r1a_virt_ensemble->get_ensemble_size(), 54);
		TS_ASSERT_EQUALS(r1a_virt_ensemble->get_conformer_table()[1]->get_nobs(), 1);
		TS_ASSERT_DELTA(r1a_virt_ensemble->get_conformer_table()[1]->get_frequency(), 0.0185185, 1.0e-6);
		TS_ASSERT( !(r1a_virt_ensemble->get_conformer_table()[1]->has_clash()) );

		NMRDummySpinlabelAtomTable & atom_table = r1a_virt_ensemble->get_conformer_table()[1]->get_atom_table();

		TS_ASSERT_DELTA(atom_table[ "N"  ].get_coordinates().x(), -0.684474, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "N"  ].get_coordinates().y(), -1.247354, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "N"  ].get_coordinates().z(), -0.321885, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "CA" ].get_coordinates().x(),  0.000000, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "CA" ].get_coordinates().y(),  0.000000, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "CA" ].get_coordinates().z(),  0.000000, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C"  ].get_coordinates().x(),  1.424894, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C"  ].get_coordinates().y(),  0.000000, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C"  ].get_coordinates().z(), -0.536021, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "O"  ].get_coordinates().x(),  2.058021, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "O"  ].get_coordinates().y(), -1.052182, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "O"  ].get_coordinates().z(), -0.646350, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "CB" ].get_coordinates().x(),  0.000000, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "CB" ].get_coordinates().y(),  0.000000, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "CB" ].get_coordinates().z(),  1.529778, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "SG" ].get_coordinates().x(),  1.036841, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "SG" ].get_coordinates().y(),  1.277633, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "SG" ].get_coordinates().z(),  2.279798, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "SD" ].get_coordinates().x(),  2.955457, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "SD" ].get_coordinates().y(),  0.548182, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "SD" ].get_coordinates().z(),  2.083698, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "CE" ].get_coordinates().x(),  3.213880, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "CE" ].get_coordinates().y(), -0.411410, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "CE" ].get_coordinates().z(),  3.614729, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "N1" ].get_coordinates().x(),  6.677743, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "N1" ].get_coordinates().y(), -1.562472, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "N1" ].get_coordinates().z(),  4.485790, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "O1" ].get_coordinates().x(),  7.741187, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "O1" ].get_coordinates().y(), -1.516146, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "O1" ].get_coordinates().z(),  5.037162, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C2" ].get_coordinates().x(),  5.645357, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C2" ].get_coordinates().y(), -0.570356, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C2" ].get_coordinates().z(),  4.681942, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C3" ].get_coordinates().x(),  4.543366, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C3" ].get_coordinates().y(), -1.088179, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C3" ].get_coordinates().z(),  3.787822, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C4" ].get_coordinates().x(),  4.912894, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C4" ].get_coordinates().y(), -2.220663, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C4" ].get_coordinates().z(),  3.174540, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C5" ].get_coordinates().x(),  6.315148, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C5" ].get_coordinates().y(), -2.619058, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C5" ].get_coordinates().z(),  3.572714, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C6" ].get_coordinates().x(),  6.316141, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C6" ].get_coordinates().y(), -3.969899, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C6" ].get_coordinates().z(),  4.284135, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C7" ].get_coordinates().x(),  7.245070, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C7" ].get_coordinates().y(), -2.621249, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C7" ].get_coordinates().z(),  2.360253, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C8" ].get_coordinates().x(),  5.198886, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C8" ].get_coordinates().y(), -0.545483, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C8" ].get_coordinates().z(),  6.141754, 1.0e-2);

		TS_ASSERT_DELTA(atom_table[ "C9" ].get_coordinates().x(),  6.129444, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C9" ].get_coordinates().y(),  0.800405, 1.0e-2);
		TS_ASSERT_DELTA(atom_table[ "C9" ].get_coordinates().z(),  4.214755, 1.0e-2);

		NMRDummySpinlabelVoxelGridCOP grid_ptr = r1a_virt_ensemble->get_voxel_grid();
		if ( grid_ptr ) {
			core::Size n_unique_sl_atoms(14);
			TS_ASSERT_EQUALS(grid_ptr->GetNumberItems(), n_unique_sl_atoms*r1a_virt_ensemble->get_ensemble_size());
			TR << "Grid dimension " << grid_ptr->GetDimension() << std::endl;
			TR << "Grid resolution " << grid_ptr->GetResolution() << std::endl;
		}
	}

	/// @brief test clash filter of ensemble conformers
	void test_quick_clash_filter() {
		using namespace core::scoring::nmr;

		TR << "Testing quick clash score calculation of NMRDummySpinlabelEnsemble." << std::endl;
		core::chemical::ResidueTypeOP restype( new core::chemical::ResidueType( *((core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard"))->get_representative_type_name3("R1A")) ) );
		NMRDummySpinlabelEnsembleOP r1a_virt_ensemble = NMRDummySpinlabelEnsembleOP( new NMRDummySpinlabelEnsemble(ensemble_conformers_database_file_, *restype) );
		TS_ASSERT_EQUALS(r1a_virt_ensemble->get_ensemble_size(), 54);
		core::Size sl_site(13);
		r1a_virt_ensemble->clash_check( cam_, sl_site, 12.0);
		auto sum_calcer = []( core::Size const previous, NMRDummySpinlabelConformerOP current)
			{ return previous + core::Size(!current->has_clash()); };
		core::Size sum = std::accumulate(r1a_virt_ensemble->get_conformer_table().begin(), r1a_virt_ensemble->get_conformer_table().end(), 0, sum_calcer);
		TS_ASSERT_EQUALS(sum, 19);
	}

	/// @brief test clash filter of ensemble conformers
	void test_elaborate_clash_filter() {
		using namespace core::scoring::nmr;
		core_init_with_additional_options("-nmr:spinlabel:elaborate_rotamer_clash_check true");
		TR << "Testing elaborate clash score calculation of NMRDummySpinlabelEnsemble." << std::endl;
		core::chemical::ResidueTypeOP restype( new core::chemical::ResidueType( *((core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard"))->get_representative_type_name3("R1A")) ) );
		NMRDummySpinlabelEnsembleOP r1a_virt_ensemble = NMRDummySpinlabelEnsembleOP( new NMRDummySpinlabelEnsemble(ensemble_conformers_database_file_, *restype) );
		TS_ASSERT_EQUALS(r1a_virt_ensemble->get_ensemble_size(), 54);
		core::Size sl_site(13);
		r1a_virt_ensemble->clash_check( cam_, sl_site, 12.0);
		auto sum_calcer = []( core::Size const previous, NMRDummySpinlabelConformerOP current)
			{ return previous + core::Size(!current->has_clash()); };
		core::Size sum = std::accumulate(r1a_virt_ensemble->get_conformer_table().begin(), r1a_virt_ensemble->get_conformer_table().end(), 0, sum_calcer);
		TS_ASSERT_EQUALS(sum, 19);

		sl_site = 12;
		r1a_virt_ensemble->clash_check( cam_, sl_site, 12.0);
		sum = std::accumulate(r1a_virt_ensemble->get_conformer_table().begin(), r1a_virt_ensemble->get_conformer_table().end(), 0, sum_calcer);
		TS_ASSERT_EQUALS(sum, 1);
		TS_ASSERT(!r1a_virt_ensemble->get_conformer_table()[1]->has_clash());
	}
};

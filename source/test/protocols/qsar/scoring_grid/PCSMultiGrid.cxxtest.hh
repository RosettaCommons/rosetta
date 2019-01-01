// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/qsar/scoring_grid/PCSMultiGrid.cxxtest.hh
/// @brief   unit test for class PCSMultiGrid which is a container for multiple PCS scoring grids
/// @details Last modified: 05/25/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/qsar/scoring_grid/PCSMultiGrid.hh>
#include <protocols/qsar/scoring_grid/PCSSingleGrid.hh>
#include <protocols/qsar/qsarMap.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <cmath>

static basic::Tracer TR("protocols.qsar.scoring_grid.PCSMultiGrid.cxxtest");

class PCSMultiGridTest : public CxxTest::TestSuite {

private:
	core::pose::Pose sh2_;

public:

	/// @brief Setup Test
	void setUp() {
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(sh2_, "protocols/qsar/scoring_grid/1x0n_model1_renum.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		sh2_.clear();
	}

	/// @brief Test instantiation of PCSMultiGrid
	void test_PCSMultiGrid_creation() {
		using namespace protocols::qsar::scoring_grid;

		// Test constructor and setup of cartesian grid
		PCSMultiGridOP grid( new PCSMultiGrid("protocols/qsar/scoring_grid/sh2_pcs_data_input.txt", 2.0));
		grid->initialize(core::Vector(-2.0, -11.0, -6.0), 10.0, 2.0);

		TS_ASSERT_EQUALS(grid->get_number_pcs_grids(), 2);
		TS_ASSERT_EQUALS(grid->get_type(), "PCSMultiGrid");
		TS_ASSERT_EQUALS(grid->grid_name(), "PCSMultiGrid");
		TS_ASSERT_DELTA(grid->get_weight(), 2.0, 1.0e-6);

		PCSSingleGridOP tb_grid = utility::pointer::dynamic_pointer_cast< PCSSingleGrid >(grid->get_pcs_grids()[1]);
		PCSSingleGridOP tm_grid = utility::pointer::dynamic_pointer_cast< PCSSingleGrid >(grid->get_pcs_grids()[2]);

		TS_ASSERT(tb_grid);
		TS_ASSERT(tm_grid);

		TS_ASSERT_DELTA(tb_grid->get_tensor()->get_ax(), -20.0, 1.0e-3);
		TS_ASSERT_DELTA(tb_grid->get_tensor()->get_rh(), -10.0, 1.0e-3);
		TS_ASSERT_DELTA(tb_grid->get_tensor()->get_alpha(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(tb_grid->get_tensor()->get_beta(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(tb_grid->get_tensor()->get_gamma(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(tb_grid->get_tensor()->get_metal_center().x(), -12.0, 1.0e-3);
		TS_ASSERT_DELTA(tb_grid->get_tensor()->get_metal_center().y(), -3.0, 1.0e-3);
		TS_ASSERT_DELTA(tb_grid->get_tensor()->get_metal_center().z(), -6.0, 1.0e-3);

		TS_ASSERT_DELTA(tm_grid->get_tensor()->get_ax(),  10.0, 1.0e-3);
		TS_ASSERT_DELTA(tm_grid->get_tensor()->get_rh(),   5.0, 1.0e-3);
		TS_ASSERT_DELTA(tm_grid->get_tensor()->get_alpha(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(tm_grid->get_tensor()->get_beta(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(tm_grid->get_tensor()->get_gamma(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(tm_grid->get_tensor()->get_metal_center().x(), -12.0, 1.0e-3);
		TS_ASSERT_DELTA(tm_grid->get_tensor()->get_metal_center().y(), -3.0, 1.0e-3);
		TS_ASSERT_DELTA(tm_grid->get_tensor()->get_metal_center().z(), -6.0, 1.0e-3);

		grid->refresh(sh2_, core::Vector(-2.0, -11.0, -6.0));
		TS_ASSERT_DELTA(tb_grid->get_pcs_values()[1]->get_pcs_exp(), 3.34975, 1.0e-3);
		TS_ASSERT_DELTA(tm_grid->get_pcs_values()[1]->get_pcs_exp(), -1.6748795, 1.0e-3);

		// Grid is centered around coordinates of K5.
		// So we test if K5 and the neighboring residues are in the grid.
		TS_ASSERT(!grid->is_in_grid(sh2_.residue(1)));
		TS_ASSERT(grid->is_in_grid(sh2_.residue(4)));
		TS_ASSERT(grid->is_in_grid(sh2_.residue(5)));
		TS_ASSERT(!grid->is_in_grid(sh2_.residue(6)));

	}

	/// @brief Test refresh of PCSMultiGrid
	void test_PCSMultiGrid_update() {
		using namespace protocols::qsar::scoring_grid;

		PCSMultiGridOP grid( new PCSMultiGrid("protocols/qsar/scoring_grid/sh2_pcs_data_input.txt"));
		grid->initialize(core::Vector(-2.0, -11.0, -6.0), 10.0, 2.0);
		grid->refresh(sh2_, core::Vector(-2.0, -11.0, -6.0));

		// Test if grid PCS values are correct
		TS_ASSERT_DELTA(grid->get_pcs_grids()[1]->get_point(-6.0,-15.0,-10.0), 2.37197662, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_pcs_grids()[1]->get_min_value(), 0.85017967, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_pcs_grids()[1]->get_max_value(), 22.13859555, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_pcs_grids()[2]->get_point(-6.0,-15.0,-10.0), -1.18598831, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_pcs_grids()[2]->get_min_value(), -11.06929778, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_pcs_grids()[2]->get_max_value(), -0.42508984, 1.0e-3);
	}

	// @brief Test score calculation of PCSMultiGrid
	void test_PCSMultiGrid_scoring() {
		using namespace protocols::qsar;
		using namespace protocols::qsar::scoring_grid;

		PCSMultiGridOP grid( new PCSMultiGrid("protocols/qsar/scoring_grid/sh2_pcs_data_input.txt", 10.0));
		grid->initialize(core::Vector(-2.0, -11.0, -6.0), 10.0, 2.0);
		grid->refresh(sh2_, core::Vector(-2.0, -11.0, -6.0));

		// PCS grid score for residue
		qsarMapCOP map;
		core::Real k5_score = grid->score(sh2_.residue(5), 9999.0, map);
		TS_ASSERT_DELTA(k5_score, 8.837, 1.0e-2);

		// PCS grid score for one atom
		core::Real k5_h_score = grid->atom_score(sh2_.residue(5), sh2_.residue(5).atom_index("H"), map);
		TS_ASSERT_DELTA(k5_h_score, 0.123, 1.0e-2);
	}
};

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/qsar/scoring_grid/PCSSingleGrid.cxxtest.hh
/// @brief   unit test for class PCSSingleGrid which is a scoring grid with PCS values measured from one lanthanide ion
/// @details Last modified: 05/24/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/qsar/scoring_grid/PCSSingleGrid.hh>
#include <protocols/qsar/qsarMap.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/conformation/Residue.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <cmath>

static basic::Tracer TR("protocols.qsar.scoring_grid.PCSSingleGrid.cxxtest");

class PCSSingleGridTest : public CxxTest::TestSuite {

private:
	core::pose::Pose sh2_;
	core::scoring::nmr::pcs::PCSTensorOP tensor_;

public:

	/// @brief Setup Test
	void setUp() {
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(sh2_, "protocols/qsar/scoring_grid/1x0n_model1_renum.pdb", core::import_pose::PDB_file);

		utility::vector1<core::Real> tensor_vals(8);
		tensor_vals[1] = 30.0;
		tensor_vals[2] = 30.0;
		tensor_vals[3] = 30.0;
		tensor_vals[4] = -12.0;
		tensor_vals[5] = -3.0;
		tensor_vals[6] = -6.0;
		tensor_vals[7] = -20.0;
		tensor_vals[8] = -10.0;
		tensor_ = core::scoring::nmr::pcs::PCSTensorOP( new core::scoring::nmr::pcs::PCSTensor );
		tensor_->set_tensor_in_pas(tensor_vals);
		tensor_->reorder_tensor();
	}

	void tearDown() {
		sh2_.clear();
		tensor_.reset();
	}

	/// @brief Test instantiation of PCSSingleGrid
	void test_PCSSingleGrid_creation() {
		using namespace protocols::qsar::scoring_grid;

		// Test the two constructors
		PCSSingleGridOP grid1( new PCSSingleGrid("protocols/qsar/scoring_grid/sim_pcs_val_tb.txt", tensor_, 2.0));
		PCSSingleGridOP grid2( new PCSSingleGrid("protocols/qsar/scoring_grid/sim_pcs_val_tb.txt", tensor_, sh2_, 2.0));

		// Test setup of cartesian grid
		grid1->initialize(core::Vector(-2.0, -11.0, -6.0), 10.0, 2.0);
		TS_ASSERT_DELTA(grid1->get_grid().getBase().x(), -7.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_grid().getBase().y(),-16.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_grid().getBase().z(),-11.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_grid().getTop().x(),   3.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_grid().getTop().y(),  -6.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_grid().getTop().z(),  -1.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_center().x(),  -2.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_center().y(), -11.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_center().z(),  -6.0, 1.0e-6);
		TS_ASSERT_EQUALS(grid1->get_dimensions().x(), 5);
		TS_ASSERT_EQUALS(grid1->get_dimensions().y(), 5);
		TS_ASSERT_EQUALS(grid1->get_dimensions().z(), 5);
		TS_ASSERT_DELTA(grid1->get_pdb_coords(0,0,0).x(), -6.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_pdb_coords(0,0,0).y(),-15.0, 1.0e-6);
		TS_ASSERT_DELTA(grid1->get_pdb_coords(0,0,0).z(),-10.0, 1.0e-6);
		TS_ASSERT_EQUALS(grid1->get_type(), "PCSSingleGrid");
		TS_ASSERT_EQUALS(grid1->grid_name(), "PCSSingleGrid");
		TS_ASSERT_DELTA(grid1->get_weight(), 2.0, 1.0e-6);
		TS_ASSERT_DELTA(grid2->get_weight(), 2.0, 1.0e-6);

		// Grid is centered around coordinates of K5.
		// So we test if K5 and the neighboring residues are in the grid.
		TS_ASSERT(!grid1->is_in_grid(sh2_.residue(1)));
		TS_ASSERT(grid1->is_in_grid(sh2_.residue(4)));
		TS_ASSERT(grid1->is_in_grid(sh2_.residue(5)));
		TS_ASSERT(!grid1->is_in_grid(sh2_.residue(6)));

		// Access the PCS data and tensor
		TS_ASSERT_DELTA(grid2->get_pcs_values()[1]->get_pcs_exp(),  3.34975, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[2]->get_pcs_exp(),  2.90150, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[3]->get_pcs_exp(),  2.83804, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[4]->get_pcs_exp(),  2.99981, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[5]->get_pcs_exp(),  3.86671, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[6]->get_pcs_exp(),  3.47057, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[7]->get_pcs_exp(),  3.21133, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[8]->get_pcs_exp(),  3.90996, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[9]->get_pcs_exp(),  2.86701, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[10]->get_pcs_exp(), 3.44099, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[11]->get_pcs_exp(), 3.73193, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_pcs_values()[12]->get_pcs_exp(), 2.79003, 1.0e-3);

		TS_ASSERT_DELTA(grid2->get_tensor()->get_ax(), -20.0, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_tensor()->get_rh(), -10.0, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_tensor()->get_alpha(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_tensor()->get_beta(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_tensor()->get_gamma(), 30.0, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_tensor()->get_metal_center().x(), -12.0, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_tensor()->get_metal_center().y(), -3.0, 1.0e-3);
		TS_ASSERT_DELTA(grid2->get_tensor()->get_metal_center().z(), -6.0, 1.0e-3);
	}

	/// @brief Test refresh of PCSSingleGrid
	void test_PCSSingleGrid_update() {
		using namespace protocols::qsar::scoring_grid;

		// Setup grid and fill it with calculated PCS values
		PCSSingleGridOP grid( new PCSSingleGrid("protocols/qsar/scoring_grid/sim_pcs_val_tb.txt", tensor_));
		grid->initialize(core::Vector(-2.0, -11.0, -6.0), 10.0, 2.0);
		grid->refresh(sh2_, core::Vector(-2.0, -11.0, -6.0));

		// Test if grid PCS values are correct
		TS_ASSERT_DELTA(grid->get_point(-6.0,-15.0,-10.0), 2.37197662, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_point(-6.0,-15.0,-2.0), 1.82975284, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_point(-2.0,-11.0,-10.0), 3.77560539, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_point(2.0,-11.0,-6.0), 1.90779777, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_point(2.0,-7.0,-2.0), 0.85017967, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_min_value(), 0.85017967, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_max_value(), 22.13859555, 1.0e-3);

	}

	// @brief Test score calculation of PCSSingleGrid
	void test_PCSSingleGrid_scoring() {
		using namespace protocols::qsar;
		using namespace protocols::qsar::scoring_grid;

		// Setup grid and fill it with calculated PCS values
		PCSSingleGridOP grid( new PCSSingleGrid("protocols/qsar/scoring_grid/sim_pcs_val_tb.txt", tensor_));
		grid->initialize(core::Vector(-2.0, -11.0, -6.0), 10.0, 2.0);
		grid->refresh(sh2_, core::Vector(-2.0, -11.0, -6.0));

		// PCS grid score for residue
		qsarMapCOP map;
		core::Real k5_score = grid->score(sh2_.residue(5), 9999.0, map);
		TS_ASSERT_DELTA(k5_score, 0.706976, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_point(-1.793,-10.890,-5.978) - grid->get_pcs_values()[5]->get_pcs_exp(),0.0993561, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_point(-0.924,-10.606,-6.332) - grid->get_pcs_values()[6]->get_pcs_exp(),-0.736973, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_point(-0.530,-11.026,-7.672) - grid->get_pcs_values()[7]->get_pcs_exp(),-0.276665, 1.0e-3);
		TS_ASSERT_DELTA(grid->get_point(-1.021,-10.032,-8.720) - grid->get_pcs_values()[8]->get_pcs_exp(), 0.278265, 1.0e-3);

		// PCS grid score for one atom
		core::Real k5_h_score = grid->atom_score(sh2_.residue(5), sh2_.residue(5).atom_index("H"), map);
		TS_ASSERT_DELTA(k5_h_score, 0.0098716346, 1.0e-3);
		TS_ASSERT_DELTA(std::sqrt(k5_h_score), grid->get_point(-1.793,-10.890,-5.978) - grid->get_pcs_values()[5]->get_pcs_exp(), 1.0e-3);
	}

};

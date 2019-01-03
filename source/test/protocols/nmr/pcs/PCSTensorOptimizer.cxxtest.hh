// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSTensorOptimizer.cxxtest.hh
/// @brief   unit test for class PCSTensorOptimizer
/// @details Last modified: 07/17/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSMultiSet.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <protocols/nmr/pcs/PCSTensorOptimizer.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/optimization/types.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>

static basic::Tracer TR("protocols.nmr.pcs.PCSTensorOptimizer.cxxtest");

class PCSTensorOptimizerTests : public CxxTest::TestSuite {
private:
	core::pose::Pose sh2_;
	utility::vector1< core::scoring::nmr::pcs::PCSSingleSetOP > pcs_data_all_lanthanides_;


public:
	/// @brief Setup Test
	void setUp() {
		using namespace core::io::nmr;
		using namespace core::scoring::nmr::pcs;

		// Initialize core & options system
		core_init();

		// Load pose from pdb
		core::import_pose::pose_from_file(sh2_, "protocols/nmr/pcs/1x0n_model1_renum.pdb", core::import_pose::PDB_file);

		// Create PCSSingleSet vector
		utility::vector1< PCSSingleSet > singleset_vec;
		PCSSingleSetOP single_dataset_tb( new PCSSingleSet( "protocols/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_ ) );
		pcs_data_all_lanthanides_.push_back(single_dataset_tb);

		PCSSingleSetOP single_dataset_dy( new PCSSingleSet( "protocols/nmr/pcs/sim_pcs_val_dy.txt", "Dy", sh2_ ) );
		pcs_data_all_lanthanides_.push_back(single_dataset_dy);

		PCSSingleSetOP single_dataset_tm( new PCSSingleSet( "protocols/nmr/pcs/sim_pcs_val_tm.txt", "Tm", sh2_ ) );
		pcs_data_all_lanthanides_.push_back(single_dataset_tm);

		PCSSingleSetOP single_dataset_er( new PCSSingleSet( "protocols/nmr/pcs/sim_pcs_val_er.txt", "Er", sh2_ ) );
		pcs_data_all_lanthanides_.push_back(single_dataset_er);
	}

	void tearDown() {
		sh2_.clear();
		pcs_data_all_lanthanides_.clear();
	}

	void test_PCSTensorOptimizer_instantiation() {
		using namespace protocols::nmr::pcs;
		PCSTensorOptimizer optimizer(pcs_data_all_lanthanides_);
	}

	void test_PCSTensorOptimizer_error_function() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;
		using namespace protocols::nmr::pcs;
		using core::pose::named_atom_id_to_atom_id;

		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("CA", 95), sh2_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("CB", 95), sh2_) );
		NMRGridSearch grid_searcher(grid_atom1, grid_atom2, sh2_, 5, 4, 0, 24);

		core::Size number_metal_ions(pcs_data_all_lanthanides_.size());

		// Set grid search center coordinates
		grid_searcher.set_grid_search_center(sh2_);

		// Set spin coordinates for every PCSSingleSet
		for ( core::Size i = 1; i <= number_metal_ions; ++i ) {
			pcs_data_all_lanthanides_[i]->update_spin_coordinates(sh2_);
		}

		// Compute combined score and PCSTensors for every PCSSingleSet
		numeric::xyzVector<core::Real> metal_coords;
		numeric::xyzVector<core::Real> best_metal_coords;
		core::Real score_svd(0);
		core::Real best_score_svd(std::numeric_limits<core::Real>::max());
		utility::vector1< core::scoring::nmr::pcs::PCSTensorCOP > singleset_tensors(number_metal_ions);

		// Grid search and SVD
		while ( grid_searcher.valid_next_grid_point(metal_coords) ) {
			score_svd = 0;
			for ( core::Size i = 1; i <= number_metal_ions; ++i ) {
				score_svd += pcs_data_all_lanthanides_[i]->solve_tensor_and_compute_score_by_svd(metal_coords) * pcs_data_all_lanthanides_[i]->get_weight();
				if ( score_svd > best_score_svd ) {
					break;
				}
			}
			if ( score_svd < best_score_svd ) {
				best_score_svd = score_svd;
				best_metal_coords = metal_coords;
				grid_searcher.set_best_grid_point(best_metal_coords);
			}
		}

		// Set PCSTensors to best values as calculated by SVD
		for ( core::Size i = 1; i <= number_metal_ions; ++i ) {
			score_svd = pcs_data_all_lanthanides_[i]->solve_tensor_and_compute_score_by_svd(best_metal_coords) * pcs_data_all_lanthanides_[i]->get_weight();
			singleset_tensors[i] = pcs_data_all_lanthanides_[i]->get_tensor_const();
		}

		TR.Debug << "Overall score calculated for group of PCSSingleSets by SVD " << best_score_svd << std::endl;
		TR.Debug << "PCS score before PCSTensorOptimizer " << score_svd << std::endl;

		core::optimization::Multivec tensor_params_all(number_metal_ions*5 + 3);
		tensor_params_all[1] = singleset_tensors[1]->get_metal_center().x();
		tensor_params_all[2] = singleset_tensors[1]->get_metal_center().y();
		tensor_params_all[3] = singleset_tensors[1]->get_metal_center().z();
		for ( core::Size i = 1; i <= singleset_tensors.size(); ++i ) {
			tensor_params_all[3 + 5 * (i-1) + 1 ] = singleset_tensors[i]->get_T_xx();
			tensor_params_all[3 + 5 * (i-1) + 2 ] = singleset_tensors[i]->get_T_xy();
			tensor_params_all[3 + 5 * (i-1) + 3 ] = singleset_tensors[i]->get_T_xz();
			tensor_params_all[3 + 5 * (i-1) + 4 ] = singleset_tensors[i]->get_T_yy();
			tensor_params_all[3 + 5 * (i-1) + 5 ] = singleset_tensors[i]->get_T_yz();
		}

		PCSTensorOptimizer optimizer(pcs_data_all_lanthanides_);
		core::Real unoptimized_score = optimizer(tensor_params_all);
		TR.Debug << "Overall score calculated by PCSTensorOptimizer " << unoptimized_score << std::endl;
		TS_ASSERT_DELTA(best_score_svd, 1.81861, 1e-2);
		TS_ASSERT_DELTA(unoptimized_score, 1.81861, 1e-2);
		TS_ASSERT_DELTA(best_score_svd,unoptimized_score, 1e-2);

		TS_ASSERT_DELTA(tensor_params_all[1], -15.862, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[2],  -3.989, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[3],  -3.266, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[4], -15.203, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[5],  -7.320, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[6],   8.681, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[7],  20.506, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[8],  12.234, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[9],  -1.018, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[10], -7.138, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[11],  7.889, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[12], 13.803, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[13], 19.657, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[14], 13.415, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[15],  3.997, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[16], -9.240, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[17],-16.339, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[18], -9.159, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[19],  2.464, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[20],  3.626, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[21], -4.737, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[22], -3.787, 1e-2);
		TS_ASSERT_DELTA(tensor_params_all[23], -3.692, 1e-2);

	}

	void test_minimization_with_PCSTensorOptimizer_() {
		using namespace protocols::nmr::pcs;
		using namespace core::scoring::nmr::pcs;

		// Create PCSTensorOptimizer
		PCSTensorOptimizer optimizer(pcs_data_all_lanthanides_);

		// Tensor values as determined in previous function
		utility::vector1< core::Real > tensors_to_optimize(pcs_data_all_lanthanides_.size()*5 + 3);
		TS_ASSERT_EQUALS(tensors_to_optimize.size(), 23);
		tensors_to_optimize[1]  = -15.862; // metal coordinates
		tensors_to_optimize[2]  =  -3.989;
		tensors_to_optimize[3]  =  -3.266;
		tensors_to_optimize[4]  = -15.203; // Tb tensor parameters
		tensors_to_optimize[5]  =  -7.320;
		tensors_to_optimize[6]  =   8.681;
		tensors_to_optimize[7]  =  20.506;
		tensors_to_optimize[8]  =  12.234;
		tensors_to_optimize[9]  =  -1.018; // Dy tensor parameters
		tensors_to_optimize[10] =  -7.138;
		tensors_to_optimize[11] =   7.889;
		tensors_to_optimize[12] =  13.803;
		tensors_to_optimize[13] =  19.657;
		tensors_to_optimize[14] =  13.415; // Tm tensor parameters
		tensors_to_optimize[15] =   3.997;
		tensors_to_optimize[16] =  -9.240;
		tensors_to_optimize[17] =  -16.339;
		tensors_to_optimize[18] =  -9.159;
		tensors_to_optimize[19] =   2.464; // Er tensor parameters
		tensors_to_optimize[20] =   3.626;
		tensors_to_optimize[21] =  -4.737;
		tensors_to_optimize[22] =  -3.787;
		tensors_to_optimize[23] =  -3.692;

		numeric::xyzVector< core::Real > metal_coord(tensors_to_optimize[1], tensors_to_optimize[2], tensors_to_optimize[3]);

		// Set spin coordinates for every PCSSingleSet in the PCSMultiSet
		for ( core::Size i = 1; i <= pcs_data_all_lanthanides_.size(); ++i ) {
			pcs_data_all_lanthanides_[i]->update_spin_coordinates(sh2_);
			pcs_data_all_lanthanides_[i]->update_matrix_A(metal_coord);
		}

		// Score before optimization
		core::Real unoptimized_score = optimizer( tensors_to_optimize );
		TR.Debug << "PCS Score before optimization " << unoptimized_score << std::endl;

		// Score and metal position after optimization
		core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.00001, true, false, false );
		core::optimization::Minimizer minimizer(optimizer, options );
		core::Real optimized_score = minimizer.run( tensors_to_optimize );
		TR.Debug << "PCS Score after optimization  " << optimized_score << std::endl;
		TR.Debug << "Optimized metal coordinates   " << tensors_to_optimize[1] << " " << tensors_to_optimize[2] << " " << tensors_to_optimize[3] << std::endl;

		TS_ASSERT_DELTA(unoptimized_score, 1.81861, 1e-2);
		TS_ASSERT_DELTA(optimized_score, 0.0000251768, 1e-2);
		TS_ASSERT_LESS_THAN(optimized_score, unoptimized_score);
	}
};

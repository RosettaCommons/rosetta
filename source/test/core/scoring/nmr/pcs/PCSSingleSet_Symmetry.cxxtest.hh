// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSSingleSet_Symmetry.cxxtest.hh
/// @brief   unit test for class PCSSingleSet that stores and handles data of one single PCS dataset (i.e. of one lanthanide ion)
///          test pcs calculation for symmetric pose either by using automatic symmetry deduction or by specifying equivalent
///          spins explicitly in the pcs input file
/// @details Last modified: 10/03/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>

static basic::Tracer TR("core.scoring.nmr.pcs.PCSSingleSet_Symmetry.cxxtest");

class PCSSingleSetSymmetryTests : public CxxTest::TestSuite {

private:
	core::pose::Pose il10_;

public:

	/// @brief Setup Test
	void setUp() {
		using namespace core::pose;
		using namespace core::conformation::symmetry;

		// Initialize core & options system
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file( il10_, "core/scoring/nmr/pcs/2ilk_chainA.pdb", core::import_pose::PDB_file);

		// Make a symmetric pose
		SymmDataOP symmdata = SymmDataOP( new SymmData );
		symmdata->read_symmetry_data_from_file( "core/scoring/nmr/pcs/2ilk_chainAB.symm" );
		core::pose::symmetry::make_symmetric_pose( il10_, *symmdata );
		runtime_assert( core::pose::symmetry::is_symmetric( il10_ ) );
	}

	void tearDown() {
		il10_.clear();
	}


	/// @brief Test PCS tensor calculation with grid search and SVD on a symmetric protein.
	///        In addition, residues in the pdb and pcs data files do not start with 1,
	///        to test that conversion from pdb to pose numbering works properly.
	///        The symmetric spins get deduced automatically from the symmetry info object
	///        and only the spins of the ASU are provided in the input file.
	void test_calc_tensor_and_score_svd_automatic_symmetry_deduction() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;
		using core::pose::named_atom_id_to_atom_id;
		using utility::to_string;

		// Important: use this option to turn automatic symmetry handling on
		core_init_with_additional_options("-nmr:pcs:use_symmetry_calc");

		core::Real weight(1.0);
		utility::vector1< PCSSingleSetOP > singleset_data;

		// Setup PCS data
		singleset_data.push_back( PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_asu_tb.txt", "Tb", il10_, weight, "CONST", "SVD" )));
		singleset_data.push_back( PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_asu_dy.txt", "Dy", il10_, weight, "CONST", "SVD" )));
		singleset_data.push_back( PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_asu_tm.txt", "Tm", il10_, weight, "CONST", "SVD" )));

		// Convert tagging residue from pdb to pose numbering
		core::pose::PDBInfoCOP pdbinfo = il10_.pdb_info();
		core::Size tag_residue(14);
		char chain_id = 'A';
		if ( pdbinfo ) {
			TR.Debug << "Converting PCS tagging residue " << to_string(tag_residue) << " " << to_string(chain_id) << " from PDB to pose numbering." << std::endl;
			tag_residue = pdbinfo->pdb2pose(chain_id, tag_residue);
		} else {
			TR.Debug << "Cannot convert PCS tagging residue from PDB to pose numbering. Assume that residue is in pose numbering instead." << std::endl;
		}

		// create gridsearch object
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("CA", tag_residue), il10_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("CB", tag_residue), il10_) );
		NMRGridSearchOP grid_searcher = NMRGridSearchOP( new NMRGridSearch(grid_atom1, grid_atom2, il10_, 5, 1, 0, 25));

		// set grid search center coordinates
		grid_searcher->set_grid_search_center(il10_);

		// Get the relevant spin coordinates from the pose
		for ( core::Size i = 1; i <= singleset_data.size(); ++i ) {
			singleset_data[i]->update_spin_coordinates(il10_);
		}

		// Run the gridsearch and score calculation
		core::Real score;
		core::Real best_score_asu(std::numeric_limits<core::Real>::max());
		numeric::xyzVector<core::Real> metal_coords;
		numeric::xyzVector<core::Real> best_metal_coords_asu;
		while ( grid_searcher->valid_next_grid_point(metal_coords) ) {
			score = 0;
			for ( core::Size i = 1; i <= singleset_data.size(); ++i ) {
				score += singleset_data[i]->solve_tensor_and_compute_score_by_svd(metal_coords) * singleset_data[i]->get_weight();
				if ( score > best_score_asu ) {
					break;
				}
			}
			if ( score < best_score_asu ) {
				best_score_asu = score;
				best_metal_coords_asu = metal_coords;
				grid_searcher->set_best_grid_point(best_metal_coords_asu);
			}
		}

		// Best metal coordinates under given grid search parameters
		TS_ASSERT_DELTA(best_metal_coords_asu.x(), 2.30373, 1e-1);
		TS_ASSERT_DELTA(best_metal_coords_asu.y(), 41.8945, 1e-1);
		TS_ASSERT_DELTA(best_metal_coords_asu.z(), 58.3661, 1e-1);

		// Show scores and determined tensors
		core::Real score1 = singleset_data[1]->solve_tensor_and_compute_score_by_svd(best_metal_coords_asu) * singleset_data[1]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[1]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score1 << std::endl;
		PCSTensorOP tensor1_asu(singleset_data[1]->get_tensor());
		tensor1_asu->show_tensor_stats(TR.Debug);
		tensor1_asu->diagonalize_tensor();
		tensor1_asu->reorder_tensor();
		TS_ASSERT_DELTA(tensor1_asu->get_T_xx(),    -17.291, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_T_xy(),     -4.969, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_T_xz(),      1.549, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_T_yy(),     17.273, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_T_yz(),      7.454, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_Eig_xx(), -2.11942, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_Eig_yy(), -18.3814, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_Eig_zz(),  20.5008, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_ax(),       30.751, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_rh(),       16.262, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_alpha(),    96.669, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_beta(),     70.572, 1e-1);
		TS_ASSERT_DELTA(tensor1_asu->get_gamma(),     9.137, 1e-1);
		tensor1_asu->show_tensor_stats(TR.Debug);

		core::Real score2 = singleset_data[2]->solve_tensor_and_compute_score_by_svd(best_metal_coords_asu) * singleset_data[2]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[2]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score2 << std::endl;
		PCSTensorOP tensor2_asu(singleset_data[2]->get_tensor());
		tensor2_asu->show_tensor_stats(TR.Debug);
		tensor2_asu->diagonalize_tensor();
		tensor2_asu->reorder_tensor();
		TS_ASSERT_DELTA(tensor2_asu->get_T_xx(),    -12.822, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_T_xy(),    -10.814, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_T_xz(),      3.210, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_T_yy(),     16.744, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_T_yz(),      5.036, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_Eig_xx(), -2.95106, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_Eig_yy(), -17.9188, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_Eig_zz(),  20.8698, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_ax(),       31.305, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_rh(),       14.968, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_alpha(),   106.981, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_beta(),     81.108, 1e-1);
		TS_ASSERT_DELTA(tensor2_asu->get_gamma(),    18.942, 1e-1);
		tensor2_asu->show_tensor_stats(TR.Debug);

		core::Real score3 = singleset_data[3]->solve_tensor_and_compute_score_by_svd(best_metal_coords_asu) * singleset_data[3]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[3]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score3 << std::endl;
		PCSTensorOP tensor3_asu(singleset_data[3]->get_tensor());
		tensor3_asu->show_tensor_stats(TR.Debug);
		tensor3_asu->diagonalize_tensor();
		tensor3_asu->reorder_tensor();
		TS_ASSERT_DELTA(tensor3_asu->get_T_xx(),     11.560, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_T_xy(),      7.608, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_T_xz(),      3.932, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_T_yy(),    -12.500, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_T_yz(),     -5.416, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_Eig_xx(),   2.7442, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_Eig_yy(),  14.1978, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_Eig_zz(),  -16.942, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_ax(),      -25.413, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_rh(),      -11.454, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_alpha(),   107.660, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_beta(),     70.438, 1e-1);
		TS_ASSERT_DELTA(tensor3_asu->get_gamma(),   168.529, 1e-1);
		tensor3_asu->show_tensor_stats(TR.Debug);

	}

	/// @brief Test PCS tensor calculation with grid search and SVD on a symmetric protein.
	///        In addition, residues in the pdb and pcs data files do not start with 1,
	///        to test that conversion from pdb to pose numbering works properly.
	///        Here, the symmetric spins are explicitly declared in the input file
	///        and the PCS is calculated by sum averaging.
	void test_calc_tensor_and_score_svd_manual_symmetry_handling() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;
		using core::pose::named_atom_id_to_atom_id;
		using utility::to_string;

		core::Real weight(1.0);
		utility::vector1< PCSSingleSetOP > singleset_data;

		// Setup PCS data
		singleset_data.push_back( PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_2su_tb.txt", "Tb", il10_, weight, "CONST", "SVD" )));
		singleset_data.push_back( PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_2su_dy.txt", "Dy", il10_, weight, "CONST", "SVD" )));
		singleset_data.push_back( PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_2su_tm.txt", "Tm", il10_, weight, "CONST", "SVD" )));
		// Turn sum averaging for the PCS of the spins from two symmetric subunits on
		for ( core::Size i = 1; i <= singleset_data.size(); ++i ) {
			singleset_data[i]->set_averaging_type("SUM");
		}

		// Convert tagging residue from pdb to pose numbering
		core::pose::PDBInfoCOP pdbinfo = il10_.pdb_info();
		core::Size tag_residue(14);
		char chain_id = 'A';
		if ( pdbinfo ) {
			TR.Debug << "Converting PCS tagging residue " << to_string(tag_residue) << " " << to_string(chain_id) << " from PDB to pose numbering." << std::endl;
			tag_residue = pdbinfo->pdb2pose(chain_id, tag_residue);
		} else {
			TR.Debug << "Cannot convert PCS tagging residue from PDB to pose numbering. Assume that residue is in pose numbering instead." << std::endl;
		}

		// create gridsearch object
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("CA", tag_residue), il10_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("CB", tag_residue), il10_) );
		NMRGridSearchOP grid_searcher = NMRGridSearchOP( new NMRGridSearch(grid_atom1, grid_atom2, il10_, 5, 1, 0, 25));

		// set grid search center coordinates
		grid_searcher->set_grid_search_center(il10_);

		// Get the relevant spin coordinates from the pose
		for ( core::Size i = 1; i <= singleset_data.size(); ++i ) {
			singleset_data[i]->update_spin_coordinates(il10_);
		}

		// Run the gridsearch and score calculation
		core::Real score;
		core::Real best_score_2su(std::numeric_limits<core::Real>::max());
		numeric::xyzVector<core::Real> metal_coords;
		numeric::xyzVector<core::Real> best_metal_coords_2su;
		while ( grid_searcher->valid_next_grid_point(metal_coords) ) {
			score = 0;
			for ( core::Size i = 1; i <= singleset_data.size(); ++i ) {
				score += singleset_data[i]->solve_tensor_and_compute_score_by_svd(metal_coords) * singleset_data[i]->get_weight();
				if ( score > best_score_2su ) {
					break;
				}
			}
			if ( score < best_score_2su ) {
				best_score_2su = score;
				best_metal_coords_2su = metal_coords;
				grid_searcher->set_best_grid_point(best_metal_coords_2su);
			}
		}

		// Best metal coordinates under given grid search parameters
		TS_ASSERT_DELTA(best_metal_coords_2su.x(), 2.30373, 1e-1);
		TS_ASSERT_DELTA(best_metal_coords_2su.y(), 41.8945, 1e-1);
		TS_ASSERT_DELTA(best_metal_coords_2su.z(), 58.3661, 1e-1);

		// Show scores and determined tensors
		core::Real score1 = singleset_data[1]->solve_tensor_and_compute_score_by_svd(best_metal_coords_2su) * singleset_data[1]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[1]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score1 << std::endl;
		PCSTensorOP tensor1_2su(singleset_data[1]->get_tensor());
		tensor1_2su->show_tensor_stats(TR.Debug);
		tensor1_2su->diagonalize_tensor();
		tensor1_2su->reorder_tensor();
		TS_ASSERT_DELTA(tensor1_2su->get_T_xx(),    -17.291, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_T_xy(),     -4.969, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_T_xz(),      1.549, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_T_yy(),     17.273, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_T_yz(),      7.454, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_Eig_xx(), -2.11942, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_Eig_yy(), -18.3814, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_Eig_zz(),  20.5008, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_ax(),       30.751, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_rh(),       16.262, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_alpha(),    96.669, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_beta(),     70.572, 1e-1);
		TS_ASSERT_DELTA(tensor1_2su->get_gamma(),     9.137, 1e-1);
		tensor1_2su->show_tensor_stats(TR.Debug);

		core::Real score2 = singleset_data[2]->solve_tensor_and_compute_score_by_svd(best_metal_coords_2su) * singleset_data[2]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[2]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score2 << std::endl;
		PCSTensorOP tensor2_2su(singleset_data[2]->get_tensor());
		tensor2_2su->show_tensor_stats(TR.Debug);
		tensor2_2su->diagonalize_tensor();
		tensor2_2su->reorder_tensor();
		TS_ASSERT_DELTA(tensor2_2su->get_T_xx(),    -12.822, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_T_xy(),    -10.814, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_T_xz(),      3.210, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_T_yy(),     16.744, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_T_yz(),      5.036, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_Eig_xx(), -2.95106, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_Eig_yy(), -17.9188, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_Eig_zz(),  20.8698, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_ax(),       31.305, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_rh(),       14.968, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_alpha(),   106.981, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_beta(),     81.108, 1e-1);
		TS_ASSERT_DELTA(tensor2_2su->get_gamma(),    18.942, 1e-1);
		tensor2_2su->show_tensor_stats(TR.Debug);

		core::Real score3 = singleset_data[3]->solve_tensor_and_compute_score_by_svd(best_metal_coords_2su) * singleset_data[3]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[3]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score3 << std::endl;
		PCSTensorOP tensor3_2su(singleset_data[3]->get_tensor());
		tensor3_2su->show_tensor_stats(TR.Debug);
		tensor3_2su->diagonalize_tensor();
		tensor3_2su->reorder_tensor();
		TS_ASSERT_DELTA(tensor3_2su->get_T_xx(),     11.560, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_T_xy(),      7.608, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_T_xz(),      3.932, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_T_yy(),    -12.500, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_T_yz(),     -5.416, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_Eig_xx(),   2.7442, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_Eig_yy(),  14.1978, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_Eig_zz(),  -16.942, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_ax(),      -25.413, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_rh(),      -11.454, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_alpha(),   107.660, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_beta(),     70.438, 1e-1);
		TS_ASSERT_DELTA(tensor3_2su->get_gamma(),   168.529, 1e-1);
		tensor3_2su->show_tensor_stats(TR.Debug);

	}

	/// @brief Test PCS tensor calculation by NLS on a symmetric protein.
	///        In addition, residues in the pdb and pcs data files do not start with 1,
	///        to test that conversion from pdb to pose numbering works properly.
	///        The symmetric spins get deduced automatically from the symmetry info object
	///        and only the spins of the ASU are provided in the input file.
	void test_calc_tensor_and_score_nls_automatic_symmetry_deduction() {
		using namespace core::scoring::nmr::pcs;
		using core::pose::named_atom_id_to_atom_id;
		using utility::to_string;

		// Important: use this option to turn automatic symmetry handling on
		// Sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pcs:use_symmetry_calc");

		core::Real weight(1.0);
		utility::vector1< PCSSingleSetOP > singleset_data;

		// Setup PCS data
		singleset_data.push_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_asu_tb.txt", "Tb", il10_, weight, "OBSIG", "NLS" )));
		singleset_data.push_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_asu_dy.txt", "Dy", il10_, weight, "OBSIG", "NLS" )));
		singleset_data.push_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_asu_tm.txt", "Tm", il10_, weight, "OBSIG", "NLS" )));

		numeric::xyzVector< core::Real > start_metal_coords(0.0, 40.0, 60.0);
		utility::fixedsizearray1<core::Real,6> metal_coords_bounds;
		metal_coords_bounds[1] = -10.0; metal_coords_bounds[2] =  30.0; metal_coords_bounds[3] =  50.0;
		metal_coords_bounds[4] =  10.0; metal_coords_bounds[5] =  50.0; metal_coords_bounds[6] =  70.0;

		for ( core::Size i = 1; i <= singleset_data.size(); ++i ) {
			singleset_data[i]->update_spin_coordinates( il10_ );
			singleset_data[i]->set_metal_coord_bounds(metal_coords_bounds);
		}

		// Score and tensor calculation
		core::Real score1 = singleset_data[1]->solve_tensor_and_compute_score_by_nls(start_metal_coords) * singleset_data[1]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[1]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score1 << std::endl;
		PCSTensorOP tensor1(singleset_data[1]->get_tensor());
		tensor1->reorder_tensor();
		tensor1->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor1->get_ax(), 27.363, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_rh(), 14.544, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_metal_center().x(),  2.327, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_metal_center().y(), 41.905, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_metal_center().z(), 57.870, 1e-1);

		core::Real score2 = singleset_data[2]->solve_tensor_and_compute_score_by_nls(start_metal_coords) * singleset_data[2]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[2]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score2 << std::endl;
		PCSTensorOP tensor2(singleset_data[2]->get_tensor());
		tensor2->reorder_tensor();
		tensor2->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor2->get_ax(), 28.283, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_rh(), 13.950, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_metal_center().x(),  2.328, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_metal_center().y(), 41.902, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_metal_center().z(), 57.873, 1e-1);

		core::Real score3 = singleset_data[3]->solve_tensor_and_compute_score_by_nls(start_metal_coords) * singleset_data[3]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[3]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score3 << std::endl;
		PCSTensorOP tensor3(singleset_data[3]->get_tensor());
		tensor3->reorder_tensor();
		tensor3->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor3->get_ax(), -22.625, 1e-1);
		TS_ASSERT_DELTA(tensor3->get_rh(), -10.014, 1e-1);
		TS_ASSERT_DELTA(tensor3->get_metal_center().x(),  2.323, 1e-1);
		TS_ASSERT_DELTA(tensor3->get_metal_center().y(), 41.906, 1e-1);
		TS_ASSERT_DELTA(tensor3->get_metal_center().z(), 57.874, 1e-1);
	}

	/// @brief Test PCS tensor calculation by NLS on a symmetric protein.
	///        In addition, residues in the pdb and pcs data files do not start with 1,
	///        to test that conversion from pdb to pose numbering works properly.
	///        Here, the symmetric spins are explicitly declared in the input file
	///        and the PCS is calculated by sum averaging.
	void test_calc_tensor_and_score_nls_manual_symmetry_handling() {
		using namespace core::scoring::nmr::pcs;
		using core::pose::named_atom_id_to_atom_id;
		using utility::to_string;

		// Set fixed RG seed for NLS fitting
		initialize_rng();

		core::Real weight(1.0);
		utility::vector1< PCSSingleSetOP > singleset_data;

		// Setup PCS data
		singleset_data.push_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_2su_tb.txt", "Tb", il10_, weight, "OBSIG", "NLS" )));
		singleset_data.push_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_2su_dy.txt", "Dy", il10_, weight, "OBSIG", "NLS" )));
		singleset_data.push_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_2ilk_2su_tm.txt", "Tm", il10_, weight, "OBSIG", "NLS" )));

		numeric::xyzVector< core::Real > start_metal_coords(0.0, 40.0, 60.0);
		utility::fixedsizearray1<core::Real,6> metal_coords_bounds;
		metal_coords_bounds[1] = -10.0; metal_coords_bounds[2] =  30.0; metal_coords_bounds[3] =  50.0;
		metal_coords_bounds[4] =  10.0; metal_coords_bounds[5] =  50.0; metal_coords_bounds[6] =  70.0;

		for ( core::Size i = 1; i <= singleset_data.size(); ++i ) {
			singleset_data[i]->update_spin_coordinates( il10_ );
			singleset_data[i]->set_metal_coord_bounds(metal_coords_bounds);
			// Turn sum averaging for the PCS of the spins from two symmetric subunits on
			singleset_data[i]->set_averaging_type("SUM");
		}

		// Score and tensor calculation
		core::Real score1 = singleset_data[1]->solve_tensor_and_compute_score_by_nls(start_metal_coords) * singleset_data[1]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[1]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score1 << std::endl;
		PCSTensorOP tensor1(singleset_data[1]->get_tensor());
		tensor1->reorder_tensor();
		tensor1->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor1->get_ax(), 27.363, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_rh(), 14.544, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_metal_center().x(),  2.327, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_metal_center().y(), 41.905, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_metal_center().z(), 57.870, 1e-1);

		core::Real score2 = singleset_data[2]->solve_tensor_and_compute_score_by_nls(start_metal_coords) * singleset_data[2]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[2]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score2 << std::endl;
		PCSTensorOP tensor2(singleset_data[2]->get_tensor());
		tensor2->reorder_tensor();
		tensor2->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor2->get_ax(), 28.283, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_rh(), 13.950, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_metal_center().x(),  2.328, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_metal_center().y(), 41.902, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_metal_center().z(), 57.873, 1e-1);

		core::Real score3 = singleset_data[3]->solve_tensor_and_compute_score_by_nls(start_metal_coords) * singleset_data[3]->get_weight();
		TR.Debug << "Calculated PCS score and tensor for dataset " << singleset_data[3]->get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score3 << std::endl;
		PCSTensorOP tensor3(singleset_data[3]->get_tensor());
		tensor3->reorder_tensor();
		tensor3->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor3->get_ax(), -22.625, 1e-1);
		TS_ASSERT_DELTA(tensor3->get_rh(), -10.014, 1e-1);
		TS_ASSERT_DELTA(tensor3->get_metal_center().x(),  2.323, 1e-1);
		TS_ASSERT_DELTA(tensor3->get_metal_center().y(), 41.906, 1e-1);
		TS_ASSERT_DELTA(tensor3->get_metal_center().z(), 57.874, 1e-1);
	}
};

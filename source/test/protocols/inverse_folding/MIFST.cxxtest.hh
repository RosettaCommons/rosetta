// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/inverse_folding/MIFST.cxxtest.hh
/// @brief  Uni tests for helper functions for using the MIF-ST model developed Yang et al.
/// @author moritzertelt (moritzertelt)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/inverse_folding/util.hh>
#include <protocols/inverse_folding/MIFST.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>


static basic::Tracer TR("MIFST");


class MIFST : public CxxTest::TestSuite {
	//Define Variables
private:
	core::pose::Pose pose_6mrr_;

public:

	void setUp() {
		core_init();
		// since only 1/3 use the pose we load it only when actually needed now
		//core::import_pose::pose_from_file(pose_6mrr_, "protocols/protein_mpnn/6MRR.pdb", core::import_pose::PDB_file);

	}

	void tearDown() {
		//pose_6mrr_.clear();
	}

	///@brief test the extract coordinates function to make sure it returns the exact coords
	void test_extract_coordinates() {
		core::import_pose::pose_from_file(pose_6mrr_, "protocols/protein_mpnn/6MRR.pdb", core::import_pose::PDB_file);
		core::select::residue_selector::ResidueSubset subset(pose_6mrr_.total_residue(), false);
		subset[1] = true;  // For this test, assuming that we are only interested in the first residue

		std::string atom_name = "CA";  // The name of the atom to extract coordinates for, here C-alpha

		std::vector<Eigen::Vector3d> coords = protocols::inverse_folding::extract_coordinates(pose_6mrr_, subset, atom_name);

		// Check the size of the returned vector
		TS_ASSERT_EQUALS(coords.size(), 1); // We only set subset[1] = true, so we should have one entry

		// Access the first residue's C-alpha coordinates from the pose for comparison
		core::conformation::Residue const & residue = pose_6mrr_.residue(1);
		core::Vector rosetta_coord = residue.xyz(atom_name);

		// Check that the coordinates match
		Eigen::Vector3d expected_coords(rosetta_coord.x(), rosetta_coord.y(), rosetta_coord.z());
		TS_ASSERT_EQUALS(coords[0], expected_coords);
	}


	/// @brief Test that calculate_Cb_position accurately calculates the Cb position
	void test_calculate_Cb_position_accuracy() {
		core::import_pose::pose_from_file(pose_6mrr_, "protocols/protein_mpnn/6MRR.pdb", core::import_pose::PDB_file);
		// Extract a sample residue from the pose (residue number 1 for example)
		core::conformation::Residue const & residue = pose_6mrr_.residue(20);

		// Extract actual atom coordinates from the residue
		core::Vector N_actual = residue.xyz("N");
		core::Vector Ca_actual = residue.xyz("CA");
		core::Vector C_actual = residue.xyz("C");
		core::Vector Cb_actual = residue.xyz("CB");

		// Convert to Eigen::Vector3d
		Eigen::Vector3d N_eigen(N_actual.x(), N_actual.y(), N_actual.z());
		Eigen::Vector3d Ca_eigen(Ca_actual.x(), Ca_actual.y(), Ca_actual.z());
		Eigen::Vector3d C_eigen(C_actual.x(), C_actual.y(), C_actual.z());

		// Calculate Cb using the original function
		Eigen::Vector3d Cb_calculated = protocols::inverse_folding::calculate_Cb_position(N_eigen, Ca_eigen, C_eigen);

		// Convert actual Cb coordinates to Eigen::Vector3d
		Eigen::Vector3d Cb_actual_eigen(Cb_actual.x(), Cb_actual.y(), Cb_actual.z());

		// Compare calculated and actual Cb positions
		core::Real tolerance = 0.1;
		TS_ASSERT_DELTA(Cb_calculated[0], Cb_actual_eigen[0], tolerance);
		TS_ASSERT_DELTA(Cb_calculated[1], Cb_actual_eigen[1], tolerance);
		TS_ASSERT_DELTA(Cb_calculated[2], Cb_actual_eigen[2], tolerance);
	}

	/// @brief Test that compute_distance_matrix accurately calculates distance matrix
	void test_compute_distance_matrix() {
		// Define some coordinates to test
		std::vector<Eigen::Vector3d> coordinates = {
			Eigen::Vector3d(0, 0, 0),
			Eigen::Vector3d(1, 0, 0),
			Eigen::Vector3d(0, 1, 0),
			Eigen::Vector3d(0, 0, 1)
			};

		// Calculate the distance matrix using the function
		Eigen::MatrixXd dists = protocols::inverse_folding::compute_distance_matrix(coordinates);

		// Verify the size of the distance matrix
		TS_ASSERT_EQUALS(dists.rows(), 4);
		TS_ASSERT_EQUALS(dists.cols(), 4);

		// Verify the diagonal elements are zero
		for ( Eigen::Index i = 0; i < 4; ++i ) {
			TS_ASSERT_DELTA(dists(i, i), 0.0, 1e-6);
		}

		// Verify other elements
		// Distance between point 0 and others
		TS_ASSERT_DELTA(dists(0, 1), 1.0, 1e-6);
		TS_ASSERT_DELTA(dists(0, 2), 1.0, 1e-6);
		TS_ASSERT_DELTA(dists(0, 3), 1.0, 1e-6);

		// Distance between point 1 and others
		TS_ASSERT_DELTA(dists(1, 2), std::sqrt(2.0), 1e-6);
		TS_ASSERT_DELTA(dists(1, 3), std::sqrt(2.0), 1e-6);

		// Distance between point 2 and others
		TS_ASSERT_DELTA(dists(2, 3), std::sqrt(2.0), 1e-6);

		// Test symmetry
		for ( Eigen::Index i = 0; i < 4; ++i ) {
			for ( Eigen::Index j = 0; j < 4; ++j ) {
				TS_ASSERT_DELTA(dists(i, j), dists(j, i), 1e-6);
			}
		}
	}

	///@brief test the get_angles function against a value computed with the original python function
	void test_get_angles() {
		// Define some coordinates for testing
		Eigen::Vector3d a1(0.0, 0.0, 0.0);
		Eigen::Vector3d b1(1.0, 0.0, 0.0);
		Eigen::Vector3d c1(1.0, 1.0, 0.0);

		Eigen::Vector3d a2(1.0, 1.0, 1.0);
		Eigen::Vector3d b2(0.0, 0.0, 0.0);
		Eigen::Vector3d c2(0.0, 1.0, 1.0);

		// Calculate the angles using the function
		core::Real calculated_angle1 = protocols::inverse_folding::get_angles(a1, b1, c1);
		core::Real calculated_angle2 = protocols::inverse_folding::get_angles(a2, b2, c2);

		// Values calculated with the original python function from yang et al.
		core::Real expected_angle1 = 1.57079633;
		core::Real expected_angle2 = 0.61547971;

		// Compare calculated value to the expected value with some tolerance
		TS_ASSERT_DELTA(calculated_angle1, expected_angle1, 1e-6);
		TS_ASSERT_DELTA(calculated_angle2, expected_angle2, 1e-6);
	}



	///@brief test the get_dihedrals function against a value computed with the original python function
	void test_get_dihedrals() {
		// Define some coordinates for testing
		Eigen::Vector3d a1(0.0, 0.0, 0.0);
		Eigen::Vector3d b1(1.0, 0.0, 0.0);
		Eigen::Vector3d c1(1.0, 1.0, 0.0);
		Eigen::Vector3d d1(0.0, 1.0, 1.0);

		Eigen::Vector3d a2(1.0, 1.0, 1.0);
		Eigen::Vector3d b2(2.0, 2.0, 2.0);
		Eigen::Vector3d c2(2.0, 1.0, 2.0);
		Eigen::Vector3d d2(1.0, 2.0, 2.0);

		// Calculate the dihedral angles using the function
		core::Real calculated_dihedral1 = protocols::inverse_folding::get_dihedrals(a1, b1, c1, d1);
		core::Real calculated_dihedral2 = protocols::inverse_folding::get_dihedrals(a2, b2, c2, d2);

		// known values from Python using the original function from yang et al.
		core::Real expected_dihedral1 = 0.78539816;
		core::Real expected_dihedral2 = -0.78539816;

		// Compare calculated values to the expected values with some tolerance
		TS_ASSERT_DELTA(calculated_dihedral1, expected_dihedral1, 1e-6);
		TS_ASSERT_DELTA(calculated_dihedral2, expected_dihedral2, 1e-6);
	}



	void test_tokenize() {
#ifndef USE_PYTORCH
		(void)protocols::inverse_folding::MASK; //avoid unused parameter in non-pytorch compile
#endif //USE_PYTORCH
		// Test case 1: Normal sequence
		std::string seq1 = "ACDE";
		core::Size mask_position1 = 1;  // Mask the second position
		std::vector<core::Size> expected1 = {0, 28, 2, 3};  // 28 is the mask character
		std::vector<core::Size> result1 = protocols::inverse_folding::tokenize(seq1, mask_position1);
		TS_ASSERT_EQUALS(result1, expected1);

		// Test case 2: Sequence with unrecognized characters
		std::string seq2 = "AC_E";
		core::Size mask_position2 = -1;  // No masking
		std::vector<core::Size> expected2 = {0, 1, 28, 3};  // 28 is mask as unrecognized should be auto masked
		std::vector<core::Size> result2 = protocols::inverse_folding::tokenize(seq2, mask_position2);
		TS_ASSERT_EQUALS(result2, expected2);
	}


	///@brief Checking some simple attributes of the process_coords function and whether it returns expected output
	void test_process_coords() {
#ifdef USE_PYTORCH
        core::import_pose::pose_from_file(pose_6mrr_, "protocols/protein_mpnn/6MRR.pdb", core::import_pose::PDB_file);
        // Create ResidueSubset according to your test
        core::select::residue_selector::ResidueSubset mock_subset(pose_6mrr_.total_residue(), true);

        // Call the process_coords function
        auto result = protocols::inverse_folding::process_coords(pose_6mrr_, mock_subset);

        torch::Tensor dist, omega, theta, phi;
        std::tie(dist, omega, theta, phi) = result;

        // Validate the dimensions of the output tensors
        TS_ASSERT_EQUALS(dist.size(0), dist.size(1));
        TS_ASSERT_EQUALS(omega.size(0), omega.size(1));
        TS_ASSERT_EQUALS(theta.size(0), theta.size(1));
        TS_ASSERT_EQUALS(phi.size(0), phi.size(1));

        // Validate that the diagonals are NaN
        for (int64_t i = 0; i < dist.size(0); ++i) {
            TS_ASSERT(std::isnan(dist[i][i].item<core::Real>()));
            TS_ASSERT(std::isnan(omega[i][i].item<core::Real>()));
            TS_ASSERT(std::isnan(theta[i][i].item<core::Real>()));
            TS_ASSERT(std::isnan(phi[i][i].item<core::Real>()));
        }

        // Validate the first few values against expected values
        torch::Tensor dist_expected = torch::tensor({
                                                            {std::nan(""), 5.654250813732024},
                                                            {5.654250813732024, std::nan("")}
                                                    });

        torch::Tensor omega_expected = torch::tensor({
                                                             {std::nan(""), 1.1788712506227315},
                                                             {1.1788712506227315, std::nan("")}
                                                     });

        torch::Tensor theta_expected = torch::tensor({
                                                             {std::nan(""), 1.4762580824529832},
                                                             {-0.5154040497802508, std::nan("")}
                                                     });

        torch::Tensor phi_expected = torch::tensor({
                                                           {std::nan(""), 0.6961051126373001},
                                                           {0.9288881336197148, std::nan("")}
                                                   });


        // checking if they match expected values, with special check for NaN since
        // it otherwise failes the TS_ASSERT even if both are NaN
        for (int64_t i = 0; i < dist_expected.size(0); ++i) {
            for (int64_t j = 0; j < dist_expected.size(1); ++j) {

                // For dist
                auto actual_dist = dist[i][j].item<core::Real>();
                auto expected_dist = dist_expected[i][j].item<core::Real>();
                if (std::isnan(actual_dist) && std::isnan(expected_dist)) continue;
                TS_ASSERT_DELTA(actual_dist, expected_dist, 1e-5);

                // For omega
                auto actual_omega = omega[i][j].item<core::Real>();
                auto expected_omega = omega_expected[i][j].item<core::Real>();
                if (std::isnan(actual_omega) && std::isnan(expected_omega)) continue;
                TS_ASSERT_DELTA(actual_omega, expected_omega, 1e-5);

                // For theta
                auto actual_theta = theta[i][j].item<core::Real>();
                auto expected_theta = theta_expected[i][j].item<core::Real>();
                if (std::isnan(actual_theta) && std::isnan(expected_theta)) continue;
                TS_ASSERT_DELTA(actual_theta, expected_theta, 1e-5);

                // For phi
                auto actual_phi = phi[i][j].item<core::Real>();
                auto expected_phi = phi_expected[i][j].item<core::Real>();
                if (std::isnan(actual_phi) && std::isnan(expected_phi)) continue;
                TS_ASSERT_DELTA(actual_phi, expected_phi, 1e-5);
            }
        }
#endif //USE_PYTORCH
	}

	void testReplaceNan() {
#ifdef USE_PYTORCH
        torch::Tensor input_tensor = torch::tensor({1.0, 2.0, std::nan("1"), 4.0});
        torch::Tensor output = protocols::inverse_folding::replace_nan(input_tensor);

        // Check if NaN is replaced
        for (auto i = 0; i < output.size(0); ++i) {
            TS_ASSERT(!std::isnan(output[i].item<core::Real>()));
        }
#endif //USE_PYTORCH
	}

	void test_prepare_sequences() {
#ifdef USE_PYTORCH
        // Test case 1: Normal sequences
        std::vector<std::string> sequences1 = {"ACDE", "FGHI"};
        std::vector<core::Size> mask_position1 = {1, 1};  // Mask the second position
        torch::Tensor expected1 = torch::tensor({{0, 28, 2, 3}, {4, 28, 6, 7}}, torch::dtype(torch::kInt64));
        torch::Tensor result1 = protocols::inverse_folding::prepare_sequences(sequences1, mask_position1);

        TS_ASSERT(torch::equal(result1, expected1));

        // Test case 2: Sequences with one unrecognized character
        std::vector<std::string> sequences2 = {"AC_E", "FGHI"};
        std::vector<core::Size> mask_position2 = {5, 5};  // No masking
        torch::Tensor expected2 = torch::tensor({{0, 1, 28, 3}, {4, 5, 6, 7}}, torch::dtype(torch::kInt64));
        torch::Tensor result2 = protocols::inverse_folding::prepare_sequences(sequences2, mask_position2);
        TS_ASSERT(torch::equal(result2, expected2));
#endif //USE_PYTORCH
	}

	void test_get_node_features() {
#ifdef USE_PYTORCH

        // Fixed example data
        torch::Tensor omega = torch::tensor({{0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}, {0.7, 0.8, 0.9}});
        torch::Tensor theta = torch::tensor({{0.2, 0.3, 0.4}, {0.5, 0.6, 0.7}, {0.8, 0.9, 1.0}});
        torch::Tensor phi = torch::tensor({{0.3, 0.4, 0.5}, {0.6, 0.7, .8}, {0.9, 1.0, 1.1}});

        // Call the C++ function
        torch::Tensor output = protocols::inverse_folding::get_node_features(omega, theta, phi);

        // Expected output from original Python function by Yang et al.
        torch::Tensor expected = torch::tensor({
                                                       {0.0000, 0.2955, 0.0000, 0.3894, 0.0000, 1.0000, 0.9553, 1.0000, 0.9211, 1.0000},
                                                       {0.1987, 0.6442, 0.4794, 0.7174, 0.5646, 0.9801, 0.7648, 0.8776, 0.6967, 0.8253},
                                                       {0.5646, 0.0000, 0.7833, 0.0000, 0.8415, 0.8253, 1.0000, 0.6216, 1.0000, 0.5403}
                                               });

        // Compare the output to the expected output
        if (!torch::allclose(output, expected, 1e-4, 1e-4)) {
            TR << "Output and expected tensors are not close!" << std::endl;
            TR << "Actual output: " << output << std::endl;
            TR << "Expected output: " << expected << std::endl;
            auto diff = torch::abs(output - expected);
            TR << "Difference: " << diff << std::endl;
            TS_ASSERT(false);
        }
#endif //USE_PYTORCH
	}
	void test_k_neighbors() {
#ifdef USE_PYTORCH
        // Create a distance tensor
        torch::Tensor dist = torch::tensor({{0.0, 1.2, 2.3, 3.4},
                                            {1.2, 0.0, 1.5, 2.6},
                                            {2.3, 1.5, 0.0, 1.7},
                                            {3.4, 2.6, 1.7, 0.0}}, torch::dtype(torch::kFloat));

        // Call the function
        torch::Tensor output = protocols::inverse_folding::get_k_neighbors(dist, 2);

        // Expected output
        torch::Tensor expected = torch::tensor({{0, 1},
                                                {1, 0},
                                                {2, 1},
                                                {3, 2}}, torch::dtype(torch::kInt64));

        // Compare the output to the expected output
        if (!torch::allclose(output, expected, 1e-4, 1e-4)) {
            TR << "Output and expected tensors are not close!" << std::endl;
            TR << "Actual output: " << output << std::endl;
            TR << "Expected output: " << expected << std::endl;
            auto diff = torch::abs(output - expected);
            TR << "Difference: " << diff << std::endl;
            TS_ASSERT(false);
        }
#endif //USE_PYTORCH
	}
	void test_get_edge_features() {
#ifdef USE_PYTORCH
        // Create example input tensors
        auto dist = torch::tensor({
                                          {0.0, 1.0, 2.0, 3.0},
                                          {1.0, 0.0, 1.0, 2.0},
                                          {2.0, 1.0, 0.0, 1.0},
                                          {3.0, 2.0, 1.0, 0.0}
                                  });

        auto omega = torch::tensor({
                                           {0.0, 0.1, 0.2, 0.3},
                                           {0.1, 0.0, 0.1, 0.2},
                                           {0.2, 0.1, 0.0, 0.1},
                                           {0.3, 0.2, 0.1, 0.0}
                                   });

        auto theta = torch::tensor({
                                           {0.0, 0.2, 0.4, 0.6},
                                           {0.2, 0.0, 0.2, 0.4},
                                           {0.4, 0.2, 0.0, 0.2},
                                           {0.6, 0.4, 0.2, 0.0}
                                   });

        auto phi = torch::tensor({
                                         {0.0, 0.3, 0.6, 0.9},
                                         {0.3, 0.0, 0.3, 0.6},
                                         {0.6, 0.3, 0.0, 0.3},
                                         {0.9, 0.6, 0.3, 0.0}
                                 });

        auto E_idx = torch::tensor({
                                           {0, 1},
                                           {1, 0},
                                           {2, 1},
                                           {3, 2}
                                   }, torch::dtype(torch::kInt64));

        // Call the C++ function
        auto output = protocols::inverse_folding::get_edge_features(dist, omega, theta, phi, E_idx);

        // Define the expected output as calculated with the original Python function from Yang et al.
        auto expected = torch::tensor({
                                              {{{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000},
                                                {1.0000, 0.0998, 0.1987, 0.1987, 0.2955, 0.2955, 0.9950, 0.9801, 0.9801, 0.9553, 0.9553}},
                                               {{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000},
                                                {1.0000, 0.0998, 0.1987, 0.1987, 0.2955, 0.2955, 0.9950, 0.9801, 0.9801, 0.9553, 0.9553}},
                                               {{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000},
                                                {1.0000, 0.0998, 0.1987, 0.1987, 0.2955, 0.2955, 0.9950, 0.9801, 0.9801, 0.9553, 0.9553}},
                                               {{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000},
                                                {1.0000, 0.0998, 0.1987, 0.1987, 0.2955, 0.2955, 0.9950, 0.9801, 0.9801, 0.9553, 0.9553}}
                                              }
                                      });

        // Compare the output to the expected output
        if (!torch::allclose(output, expected, 1e-4, 1e-4)) {
            TR << "Output and expected tensors are not close!" << std::endl;
            TR << "Actual output: " << output << std::endl;
            TR << "Expected output: " << expected << std::endl;
            auto diff = torch::abs(output - expected);
            TR << "Difference: " << diff << std::endl;
            TS_ASSERT(false);
        }
#endif //USE_PYTORCH
	}

	void test_get_mask() {
#ifdef USE_PYTORCH
        // Create example input tensor (3 nodes, 2 neighbors, 6 edge features)
        auto E = torch::tensor({
                                       {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, {7.0, 8.0, 9.0, 10.0, 11.0, 12.0}},
                                       {{1.0, 2.0, 3.0, 4.0, 5.0, std::nan("")}, {7.0, 8.0, 9.0, 10.0, 11.0, 12.0}},
                                       {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, {7.0, 8.0, 9.0, 10.0, 11.0, 12.0}}
                               }, torch::dtype( torch::kFloat ));
        // Call the C++ function
        auto output = protocols::inverse_folding::get_mask(E);
        // Define the expected output
        auto expected = torch::tensor({
                                              {{1.0}, {1.0}},
                                              {{0.0}, {1.0}},
                                              {{1.0}, {1.0}}
                                      }, torch::dtype(torch::kFloat));

        // Compare the output to the expected output
        if (!torch::allclose(output, expected, 1e-4, 1e-4)) {
            TR << "Output and expected tensors are not close!" << std::endl;
            TR << "Actual output: " << output << std::endl;
            TR << "Expected output: " << expected << std::endl;
            auto diff = torch::abs(output - expected);
            TR << "Difference: " << diff << std::endl;
            TS_ASSERT(false);
        }
#endif //USE_PYTORCH
	}

	void test_collate_structure() {
#ifdef USE_PYTORCH

        // Load the pose from the PDB file
        core::import_pose::pose_from_file(pose_6mrr_, "protocols/protein_mpnn/6MRR.pdb", core::import_pose::PDB_file);

        core::select::residue_selector::ResidueSubset subset(pose_6mrr_.total_residue(), true);

        auto result = protocols::inverse_folding::process_coords(pose_6mrr_, subset);

        torch::Tensor dist, omega, theta, phi;
        std::tie(dist, omega, theta, phi) = result;
        // Prepare data for collate_structure
        std::vector<std::string> sequences = {pose_6mrr_.sequence()};
        std::vector<core::Size> mask_positions = {500}; // no masking

        auto collated_data_cpp = protocols::inverse_folding::collate_structure( sequences,
                                                                      dist,
                                                                      omega,
                                                                      theta,
                                                                      phi,
                                                                      mask_positions);

        torch::Tensor src_cpp, nodes_cpp, edges_cpp, connections_cpp, edge_mask_cpp;
        std::tie( src_cpp, nodes_cpp, edges_cpp, connections_cpp, edge_mask_cpp ) = collated_data_cpp;


        // Define expected values from Python
        auto src_expected = torch::tensor({5, 18}, torch::dtype(torch::kInt64));
        auto nodes_expected = torch::tensor({
                                                    {{0.0, 0.9955},
                                                     {0.9242, 0.4792}}
                                            }, torch::dtype(torch::kFloat32));
        auto edges_expected = torch::tensor({{{5.6543}}}, torch::dtype(torch::kFloat32));
        auto connections_expected = torch::tensor({{{1, 30, 29, 2, 28},
                                                    {0, 2, 6, 3, 28},
                                                    {5, 4, 28, 3, 6},
                                                    {4, 2, 6, 7, 1},
                                                    {3, 5, 2, 7, 8}}}, torch::dtype(torch::kInt64));
        auto edge_mask_expected = torch::tensor({{{1.0}, {1.0}},
                                                 {{1.0}, {1.0}}}, torch::dtype(torch::kFloat32));


        // Check if the first few values of the C++ and Python results are equivalent
        TS_ASSERT(torch::allclose(src_cpp.slice(1, 0, 2), src_expected, 1e-4));
        TS_ASSERT(torch::allclose(nodes_cpp.slice(1, 0, 2).slice(2, 0, 2), nodes_expected, 1e-4));
        TS_ASSERT(torch::allclose(edges_cpp.slice(1, 0, 1).slice(2, 0, 1).slice(3, 0, 1), edges_expected, 1e-4));
        TS_ASSERT(torch::allclose(connections_cpp.slice(1, 0, 5).slice(2, 0, 5), connections_expected, 1e-4));
        TS_ASSERT(torch::allclose(edge_mask_cpp.slice(1, 0, 2).slice(2, 0, 2).slice(3, 0, 1), edge_mask_expected, 1e-4));
#endif //USE_PYTORCH
	}

	void test_raw_predict() {
#ifdef USE_PYTORCH
        // Load the pose from the PDB file
        core::import_pose::pose_from_file(pose_6mrr_, "protocols/protein_mpnn/6MRR.pdb", core::import_pose::PDB_file);

        // just predicting for the two last residues
        core::select::residue_selector::ResidueSubset residue_selection(pose_6mrr_.total_residue(), false);
        residue_selection[67] = true;
        residue_selection[66] = true;

        // while taking the whole pose into account for features
        core::select::residue_selector::ResidueSubset feature_selection(pose_6mrr_.total_residue(), true);

        torch::Tensor logits;
        // Forward pass
        logits = protocols::inverse_folding::MIFST::get_instance()->predict(pose_6mrr_, residue_selection, feature_selection, false, false);
        // Extract the relevant part of the tensor
        logits = logits.cpu();
        auto relevant_output = logits.index({0, 1}).to(torch::kCPU, torch::kFloat64).contiguous().view({-1});

        // Convert the tensor to a std::vector
        std::vector<core::Real> predicted_output(relevant_output.data_ptr<core::Real>(), relevant_output.data_ptr<core::Real>() + relevant_output.numel());

        // Expected output from Python
        std::vector<core::Real> expected_output = {
                1.2110, -1.3151, 0.5495, 1.7907, 0.7511, -0.3028, 2.8298,
                3.0496, 1.9913, 1.5270, 3.7578, 1.9784, -5.1468, 0.6510,
                2.6342, 2.7779, 1.5619, 2.4155, 2.2324, -0.1481, -14.5878,
                -20.7614, -8.7002, -21.0255, -14.5939, -11.7598, -20.7413,
                -21.4430, -19.8671, -19.9772
        };

        // Compare expected and predicted output
        for (core::Size i = 0; i < expected_output.size(); ++i) {
            TS_ASSERT_DELTA(predicted_output[i], expected_output[i], 1);
        }
#endif //USE_PYTORCH
	}

};

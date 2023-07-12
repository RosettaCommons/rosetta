// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /esm_perplexity/EsmPerplexity.cxxtest.hh
/// @brief Unit tests for the ESM language model protocol
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
# include <protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>

// std headers
#include <numeric> // needed to not break windows build

using namespace protocols;
using namespace protocols::esm_perplexity;

class EsmPerplexityTest : public CxxTest::TestSuite {

private:
	std::string existing_model_ = "esm2_t6_8M_UR50D";
	std::string wrong_model_ = "thisShouldntWork";

#ifdef USE_TENSORFLOW
    // std::shared_ptr<EsmPerplexityTensorflowProtocol> tensorflow_protocol_;
#endif // USE_TENSORFLOW

	void setUp() override {
		core_init();
#ifdef USE_TENSORFLOW
#endif // USE_TENSORFLOW
	}

	void tearDown() override {
#ifdef USE_TENSORFLOW
#endif // USE_TENSORFLOW
	}

public:

	/// @brief Check that we fail on wrong name input
	void test_protocol_fails_with_wrong_name() {
#ifdef USE_TENSORFLOW
        TS_ASSERT_THROWS_ANYTHING( EsmPerplexityTensorflowProtocol failing_protocol( wrong_model_ ) );
#endif // USE_TENSORFLOW
	}

	/// @brief make sure that the softmax function returns correct values
	/// @note we don't need USE_TENSORFLOW here since its always compiled
	void test_softmax() {

		std::map<core::Size, utility::vector1<core::Real>> logit_map = {
			{1, {1.0, 2.0, 3.0, 6.0}},
			{2, {4.0, 2.0, 9.0, 1.0}},
			{3, {-1.0, -2.0, -3.0, -4.0}},
			{4, {4e150, 2e150, 9e150, 1e150}}, // test for large numbers
			{5, {-1e-150, -2e-150, -3e-150, -4e-150}} // test for small numbers
			};

		std::map<core::Size, utility::vector1<core::Real>> expected_softmax_map = {
			{1, {0.00626879, 0.01704033, 0.04632042, 0.93037047}},
			{2, {0.00668457, 0.00090466, 0.99207797, 0.00033281}},
			{3, {0.64391426, 0.23688282, 0.08714432, 0.0320586}},
			{4, {0.0, 0.0, 1.0, 0.0}},
			{5, {0.25, 0.25, 0.25, 0.25}}
			};

		std::map<core::Size, utility::vector1<core::Real>> softmax_map;

		EsmPerplexityTensorflowProtocol::softmax(logit_map, softmax_map);

		TS_ASSERT_EQUALS(logit_map.size(), softmax_map.size()); // ensure size equality between input and output
		for ( const auto& pair : softmax_map ) {
			double sum = std::accumulate(pair.second.begin(), pair.second.end(), 0.0);
			TS_ASSERT_DELTA(sum, 1.0, 1e-6); // The sum of softmax values should be approximately 1

			for ( size_t i = 0; i < pair.second.size(); ++i ) {
				TS_ASSERT_DELTA(pair.second[i+1], expected_softmax_map[pair.first][i+1], 1e-6); // check that we match expected output
			}
		}
	}

	/// @brief make sure the feature creation (tokenization) is working correctly (even when including NCAAs)
	void test_copy_feat_to_tensor() {
#ifdef USE_TENSORFLOW
        // Construct a pose.
        core::pose::Pose pose;
        core::pose::make_pose_from_sequence( pose, "LAGVSERTIDPKQNFYMHWCX[ORN]", "fa_standard"); // AA's in the order of the token alphabet

        // Define the position, attention mask, and input tensor containers.
        core::Size position = 2;
        utility::vector1<core::Size> attention_mask_vec = {1, 3, 4};\

        // get chain of selected residue
        core::Size const chain = pose.residue(position).chain();
        // get chain sequence
        std::string const chain_sequence = pose.chain_sequence(chain);
        // get chain length
        long const chain_length = static_cast<long>(chain_sequence.size());

        // prepare input tensorflow containers
        basic::tensorflow_manager::RosettaTensorflowTensorContainer<int> input_tensor_seq;
        basic::tensorflow_manager::RosettaTensorflowTensorContainer<int> input_tensor_mask;
        input_tensor_seq.initialize(TF_INT32, {1, chain_length + 2});
        input_tensor_mask.initialize(TF_INT32, {1, chain_length + 2});
        // run function
        EsmPerplexityTensorflowProtocol::copy_feat_to_tensor(
                pose, position, attention_mask_vec, input_tensor_seq, input_tensor_mask
        );
        // should start with 0 ("cls" token), then counting up from 4 to 23 for the amino acids except for the position which is being predicted (32)
        // and the NCAA which is set to the unknown token (3)
        std::vector<int> expected_input_tensor_seq = {0, 4, 32, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 3, 2};
        // zeros should be at the positions specified by the attention_mask_vec (1-indexed)
        std::vector<int> expected_input_tensor_mask = {0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        // compare expected results with ones returned
        for (core::Size i = 0; i < expected_input_tensor_seq.size(); ++i) {
            TS_ASSERT_EQUALS(input_tensor_seq(i + 1), expected_input_tensor_seq[i]);
        }
        for (core::Size i = 0; i < expected_input_tensor_mask.size(); ++i) {
            TS_ASSERT_EQUALS(input_tensor_mask(i + 1), expected_input_tensor_mask[i]);
        }
#endif // USE_TENSORFLOW
	}

};


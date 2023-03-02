// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protein_mpnn/ProteinMPNN.cxxtest.hh
/// @brief Unit tests for ProteinMPNN sequence predictor
/// @author Frederick Chan (fredchan@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
# include <protocols/protein_mpnn/ProteinMPNN.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

class ProteinMPNNTest : public CxxTest::TestSuite {
private:
	core::pose::Pose pose_6mrr_; // 6MRR is the PDB I have exported inputs/outputs for in Python
	core::pose::Pose pose_5l33_;

	/// @brief Get sequence recovery between two sequences
	core::Real sequence_recovery( std::string const & original_seq, std::string const & sample_seq ) {
		CxxTest::setAbortTestOnFail(true); // This will reset after the test that uses this function finishes
		TS_ASSERT_EQUALS(original_seq.length(), sample_seq.length()); // Prevent undefined behavior

		int recovered = 0;
		for ( core::Size residue = 0; residue < sample_seq.size(); ++residue ) {
			if ( sample_seq[residue] == original_seq[residue] ) {
				recovered++;
			}
		}
		return ( recovered / core::Real( original_seq.size() ) );
	}

#ifdef USE_PYTORCH
	/// @brief Set options for testing deterministic version ProteinMPNN
	/// @note in case ProteinMPNNOptions defaults change, we still want these options to remain the same
	void set_deterministic_options( protocols::protein_mpnn::ProteinMPNNOptions & options, core::pose::Pose & pose ) {
		options.sequence = pose.sequence();

		options.coord_mask.clear();
		options.coord_mask.resize( pose.total_residue(), true );

		options.chain_mask.clear();
		options.chain_mask.resize( pose.num_chains(), true );

		options.pos_mask.clear();
		options.pos_mask.resize( pose.total_residue(), true );

		options.omit_AAs.clear();
		options.omit_AAs.push_back( 'X' );

		options.bias_AAs.clear();
		options.bias_AAs.resize( 21, 0.0 );

		options.temperature = 0.1;
		options.tied_positions.clear();
		options.batch_size = 1;
		options.deterministic_flag = true;
	}
#endif //USE_PYTORCH

public:
	/// @brief Setup Test
	void setUp() {
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(pose_6mrr_, "protocols/protein_mpnn/6MRR.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file(pose_5l33_, "protocols/protein_mpnn/5L33.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		pose_6mrr_.clear();
		pose_5l33_.clear();
	}

	/// @brief Test that we can initialize a ProteinMPNNOptions with a pose
	void test_ProteinMPNNOptions_init() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );

		TS_ASSERT_EQUALS( options_6mrr.sequence, pose_6mrr_.sequence() );
#endif //USE_PYTORCH
	}

	/// @brief Test that we are getting sane sequence recovery rates for sample without custom options
	/// @note Around 50% is typical; higher is better
	void test_ProteinMPNN_sample_seq_recovery_at_least_40_pct() {
#ifdef USE_PYTORCH
		std::string original_seq = pose_6mrr_.sequence();
		std::string sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_ );

		TS_ASSERT( sequence_recovery( original_seq, sample_seq ) >= 0.4 );
#endif //USE_PYTORCH
	}

	/// @brief Test that we are getting sane sequence recovery rates for sample with custom options
	/// @note Around 50% is typical; higher is better
	void test_ProteinMPNN_sample_with_opts_seq_recovery_at_least_40_pct() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );

		std::string original_seq = pose_6mrr_.sequence();
		std::string sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr );

		TS_ASSERT( sequence_recovery( original_seq, sample_seq ) >= 0.4 );
#endif //USE_PYTORCH
	}

	/// @brief Test that we are getting sane sequence recovery rates for batch sample of batch_size=2
	/// @note Around 50% is typical; higher is better
	//void test_ProteinMPNN_sample_batch_2_seq_recovery_at_least_40_pct() {
	// protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
	// options_6mrr.batch_size = 2;

	// std::string original_seq = pose_6mrr_.sequence();
	// utility::vector1< std::string > sample_seqs = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample_batch( pose_6mrr_, options_6mrr );

	// TS_ASSERT_EQUALS( sample_seqs.size(), 2 );
	// for ( auto & sample_seq : sample_seqs ) {
	//  TS_ASSERT( sequence_recovery( original_seq, sample_seq ) >= 0.4 );
	// }
	//}

	/// @brief Test that the input and output sequences are the same shapes for sample without custom options
	void test_ProteinMPNN_sample_shape_unchanged() {
#ifdef USE_PYTORCH
		std::string original_seq = pose_6mrr_.sequence();
		std::string sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_ );

		TS_ASSERT_EQUALS(original_seq.length(), sample_seq.length())
#endif //USE_PYTORCH
	}

	/// @brief Test that the input and output sequences are the same shapes for sample with custom options
	void test_ProteinMPNN_sample_with_opts_shape_unchanged() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );

		std::string original_seq = pose_6mrr_.sequence();
		std::string sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr );

		TS_ASSERT_EQUALS(original_seq.length(), sample_seq.length())
#endif //USE_PYTORCH
	}

	/// @brief Test that the input and output sequences are the same shapes for batch sample of batch_size=2
	//void test_ProteinMPNN_sample_batch_2_shape_unchanged() {
	// protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
	// options_6mrr.batch_size = 2;

	// std::string original_seq = pose_6mrr_.sequence();
	// utility::vector1< std::string > sample_seqs = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample_batch( pose_6mrr_, options_6mrr );

	// TS_ASSERT_EQUALS( sample_seqs.size(), 2 );
	// for ( auto & sample_seq : sample_seqs ) {
	//  TS_ASSERT_EQUALS( original_seq.length(), sample_seq.length() );
	// }
	//}

	/// @brief Test that setting deterministic mode results in the same sequence every time
	void test_ProteinMPNN_sample_deterministic_always_same() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		options_6mrr.deterministic_flag = true;

		std::string sample_seq_prev = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr );
		std::string sample_seq;
		for ( int i = 0; i < 4; ++i ) {
			sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr );
			TS_ASSERT_EQUALS( sample_seq, sample_seq_prev );
			sample_seq_prev = sample_seq;
		}
#endif //USE_PYTORCH
	}

	/// @brief Test that setting deterministic mode results in the same sequences as the Python model
	void test_ProteinMPNN_sample_deterministic_matches_python() {
#ifdef USE_PYTORCH
		std::string sample_seq;

		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		set_deterministic_options( options_6mrr, pose_6mrr_ );
		sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr );
		TS_ASSERT_EQUALS( sample_seq, "GVSEELEKYRKELEAFLKKEGITNVKIEIKDGELKIETKGGSEKLKKKLEELKKKLEEKGYKVEVKIE" );

		protocols::protein_mpnn::ProteinMPNNOptions options_5l33( pose_5l33_ );
		set_deterministic_options( options_5l33, pose_5l33_ );
		sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_5l33_, options_5l33 );
		TS_ASSERT_EQUALS( sample_seq, "SVPADEAKALEFVKALEEGNPELMKSVISPDTKMEVNGKEYVGEEIVEYVKEIKKKGTKYKLISYKKVGDKYLFEVEVENNGKKYKAKIEIEVKDGKIAKVIITCK" );
#endif //USE_PYTORCH
	}

	/// @brief Test that setting an invalid sequence option throws an exception
	void test_ProteinMPNN_sample_sequence_invalid_sequence_option_throws_exception() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		options_6mrr.sequence = "X";
		TS_ASSERT_THROWS_ANYTHING( protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr ) );
#endif //USE_PYTORCH
	}

	/// @brief Test that setting an invalid coord mask option throws an exception
	void test_ProteinMPNN_sample_sequence_invalid_coord_mask_option_throws_exception() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		options_6mrr.coord_mask = utility::vector1_bool();
		TS_ASSERT_THROWS_ANYTHING( protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr ) );
#endif //USE_PYTORCH
	}

	/// @brief Test that setting an invalid chain mask option throws an exception
	void test_ProteinMPNN_sample_sequence_invalid_chain_mask_option_throws_exception() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		options_6mrr.chain_mask = utility::vector1_bool();
		TS_ASSERT_THROWS_ANYTHING( protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr ) );
#endif //USE_PYTORCH
	}

	/// @brief Test that setting an invalid pos mask option throws an exception
	void test_ProteinMPNN_sample_sequence_invalid_pos_mask_option_throws_exception() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		options_6mrr.pos_mask = utility::vector1_bool();
		TS_ASSERT_THROWS_ANYTHING( protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr ) );
#endif //USE_PYTORCH
	}

	/// @brief Test that setting an invalid bias AAs option throws an exception
	void test_ProteinMPNN_sample_sequence_invalid_bias_AAs_option_throws_exception() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		options_6mrr.bias_AAs = utility::vector1< core::Real >();
		TS_ASSERT_THROWS_ANYTHING( protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr ) );
#endif //USE_PYTORCH
	}

	/// @brief Test that setting an invalid tied positions option throws an exception
	void test_ProteinMPNN_sample_sequence_invalid_tied_positions_option_throws_exception() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );

		utility::vector1< utility::vector1< core::Size > > tied_positions;
		utility::vector1< core::Size > tied_position;
		tied_position.push_back( 1 );
		tied_position.push_back( 9999 );
		tied_position.push_back( 9999 );
		tied_positions.push_back( tied_position );

		options_6mrr.tied_positions = tied_positions;
		TS_ASSERT_THROWS_ANYTHING( protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr ) );
#endif //USE_PYTORCH
	}

	/// @brief Test that pos mask option fixes the residue type at the selected residue numbers.
	void test_ProteinMPNN_sample_sequence_pos_mask_fixes_residues() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		utility::vector1< bool > pos_mask;
		pos_mask.resize( pose_6mrr_.total_residue(), true );
		for ( int residue = 1; residue <= 10; ++residue ){
			pos_mask[residue] = false;
		}
		options_6mrr.pos_mask = pos_mask;
		std::string sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr );

		std::string first_residues_sample_seq = sample_seq.substr(0, 10);
		std::string first_residues_6mrr = pose_6mrr_.sequence().substr(0, 10);

		TS_ASSERT_EQUALS(first_residues_6mrr, first_residues_sample_seq);
#endif //USE_PYTORCH
	}

	/// @brief Test that omit_AA_pos option blocks specified residue types at the selected residue numbers.
	void test_ProteinMPNN_sample_sequence_omit_AA_pos_restricts_types() {
#ifdef USE_PYTORCH
		protocols::protein_mpnn::ProteinMPNNOptions options_6mrr( pose_6mrr_ );
		utility::vector1< char > omit_except_W{ 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','X','Y' };
		for ( int residue = 1; residue <= 10; ++residue ){
			options_6mrr.omit_AAs_pos[residue] = omit_except_W;
		}
		std::string sample_seq = protocols::protein_mpnn::ProteinMPNN::get_instance()->sample( pose_6mrr_, options_6mrr );

		std::string first_residues_sample_seq = sample_seq.substr(0, 10);
		std::string first_residues_W = "WWWWWWWWWW";

		TS_ASSERT_EQUALS(first_residues_W, first_residues_sample_seq);
#endif //USE_PYTORCH
	}
};

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_mpnn/ProteinMPNN.hh
/// @brief For the ProteinMPNN PyTorch model developed by Dauparas et al.
/// @author Brian Koepnick (koepnick@uw.edu)

#ifdef USE_TORCH

#ifndef INCLUDED_protocols_protein_mpnn_ProteinMPNN_hh
#define INCLUDED_protocols_protein_mpnn_ProteinMPNN_hh

#include <protocols/protein_mpnn/ProteinMPNN.hh>
#include <basic/citation_manager/CitationCollectionBase.fwd.hh>

// Core headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/SingletonBase.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility>

#include <torch/script.h>

namespace protocols {
namespace protein_mpnn {

/// @brief Options to configure ProteinMPNN sampling
struct ProteinMPNNOptions {

	/// @brief Constructor
	/// @details Constructor requires pose to initialize options with the correct shape
	ProteinMPNNOptions( const core::pose::Pose & pose );

	/// @brief Reset options to defaults
	void
	reset( const core::pose::Pose & pose );

	// the user can override the starting sequence, why not
	std::string sequence;

	// list of bools for omitted coordinates
	// length equal to pose.total_residue()
	utility::vector1_bool coord_mask;

	// fixed positions by chain ID
	// NB: chain_mask uses chain_ID (NOT PDB chain lettering)
	utility::vector1_bool chain_mask;

	// fixed positions by residue number
	// NB: this supersedes "chain_mask_pos" from the python ProteinMPNN,
	// which uses chain-specific residue numbering
	// in Rosetta, it's usually easier to work with full-pose residue numbering
	utility::vector1_bool pos_mask;

	// global disallowed AA types (plus unknown 'X')
	utility::vector1< char > omit_AAs;

	// position-specific disallowed AA types (plus unknown 'X')
	utility::vector1< utility::vector1< char > > omit_AAs_pos;

	// global biases for AA types (plus unknown 'X')
	utility::vector1< core::Real > bias_AAs;

	// temperature for MCMC sampling
	core::Real temperature;

	// collection of lists, where each list contains tied positions, e.g.
	// NB: unlike the python ProteinMPNN, this uses full-pose residue numbering,
	// so chain ID is omitted from tied positions
	utility::vector1< utility::vector1< core::Size > > tied_positions;

	// batch_size
	// currently, only supports batch_size = 1
	core::Size batch_size;

	// whether we should sample with deterministic mode
	bool deterministic_flag;

	// TODO: implement PSSM options??
	// pssm
	// pssm_multi // pssm weight(?)
	// pssm_threshold
	// tied_positions

private:

	// default constructor disabled
	// object MUST be constructed with parameters
	ProteinMPNNOptions();
};


/// @brief Interface for ProteinMPNN neural network
/// @details Loads the pre-trained PyTorch model from a TorchScript file, and provides access to sample() method for predicting the sequence of an input protein backbone.
class ProteinMPNN : public utility::SingletonBase< ProteinMPNN >{
	friend class utility::SingletonBase< ProteinMPNN >;

	typedef torch::List< torch::List< int64_t > > TiedPositions;

private:

	/// @brief Private constructor for singleton class
	ProteinMPNN();
	ProteinMPNN( ProteinMPNN const & ) = delete;
	ProteinMPNN operator=( ProteinMPNN const & ) = delete;


public:

	/// @brief Call the ProteinMPNN sample() method
	/// @details Predict batch of sequences for an input protein backbone
	utility::vector1< std::string >
	sample_batch( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options );

	/// @brief Call the ProteinMPNN sample() method
	/// @details Predict a single sequence for an input protein backbone
	std::string
	sample( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options );

	/// @brief Call the ProteinMPNN sample() method
	/// @details Predict a single sequence for an input protein backbone, using default options
	std::string
	sample( const core::pose::Pose & pose );

    ///@brief Predict amino acid probabilities/logits and fill them into maps
    ///@details Runs ProteinMPNN inference using the provided options and fills the predicted probabilities/logits into provided maps
    ///@param[in] pose The pose which will be used for prediction
    ///@param[in] mpnn_options An ProteinMPNNOptions object specifying options for inference
    ///@param[in] probs_map A map of residue positions to vectors of amino acid probabilities that will be filled
    ///@param[in] logits_map A map of residue positions to vectors of amino acid logits that will be filled
    void
    get_probabilities_and_logits( core::pose::Pose const & pose, ProteinMPNNOptions const & mpnn_options,
                                               std::map<core::Size, utility::vector1<core::Real>>& probs_map,
                                               std::map<core::Size, utility::vector1<core::Real>>& logits_map);

	/// @brief Get the citation for ProteinMPNN
	/// @details TODO: fill in details for Dauparas et al.
	static
	basic::citation_manager::CitationCollectionBaseCOP
	get_ProteinMPNN_neural_net_citation();

private:

	/// @brief Initialize the PyTorch model used by this class.
	/// @details Called ONCE by class constructor.  Triggers read from disk!
	void init_mpnn_model();

	/// @brief Generate PyTorch tensor inputs from user options
	std::vector< torch::jit::IValue >
	generate_sample_inputs( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options, bool tied = false );

	/// @brief Create PyTorch optional inputs from config options
	std::unordered_map< std::string, torch::jit::IValue >
	generate_sample_optional_inputs( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options, bool tied = false );

	/// @brief Check option validity
	bool
	verify_options( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options );

	torch::jit::script::Module mpnn_module;

};

} //protein_mpnn
} //protocols

#endif //INCLUDED_protocols_protein_mpnn_ProteinMPNN_hh

#endif //USE_TORCH

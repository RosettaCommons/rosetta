// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_mpnn/ProteinMPNN.cc
/// @brief For the ProteinMPNN PyTorch model developed by Dauparas et al.
/// @author Brian Koepnick (koepnick@uw.edu)

#ifdef USE_TORCH

// Project headers:
#include <protocols/protein_mpnn/ProteinMPNN.hh>
#include <protocols/protein_mpnn/util.hh>

// Core headers:
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

// Basic headers:
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>

// Utility headers:
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/memory.hh>

#include <torch/csrc/autograd/generated/variable_factories.h>
#include <torch/script.h>
#include <algorithm>

static basic::Tracer TR( "protocols.protein_mpnn.ProteinMPNN" );

namespace protocols {
namespace protein_mpnn {

using namespace core;
using namespace core::pose;

/// @brief Constructor.
ProteinMPNN::ProteinMPNN()
{
	init_mpnn_model();
}

////////////////////////////////
// PUBLIC MEMBER FUNCTIONS:
////////////////////////////////

/// @brief Sample sequences using provided options
utility::vector1< std::string >
ProteinMPNN::sample_batch( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options ){
	TR.Debug << "Invoking ProteinMPNN.sample()... " << std::endl;
	utility::vector1< std::string > sample_seqs;

	// verify options
	if( !verify_options( pose, mpnn_options ) ){
		TR.Error << "Bad inputs, aborting ProteinMPNN::sample()..." << std::endl;
		throw CREATE_EXCEPTION( utility::excn::Exception, "Invalid options passed to ProteinMPNN. Check log for details." );
	}

	// do we want the model's sample() or tied_sample() method
	bool tied( false );
	torch::jit::script::Method sample_method = mpnn_module.get_method( "sample" );
	if( !mpnn_options.tied_positions.empty() ){
		tied = true;
		sample_method = mpnn_module.get_method( "tied_sample" );
	}

	// prepare PyTorch inputs
	std::vector<torch::jit::IValue> inputs = generate_sample_inputs( pose, mpnn_options, tied );
	std::unordered_map< std::string, torch::jit::IValue > optional_inputs = generate_sample_optional_inputs( pose, mpnn_options, tied );
	TR.Debug << "ProteinMPNN.sample() inputs: " << std::endl;
	for( auto & input : inputs ){
		TR.Debug << input << std::endl;
	}
	TR.Debug << "ProteinMPNN.sample() optional inputs: " << std::endl;
	for ( auto & item : optional_inputs ) {
		TR.Debug << item.first << ": " << item.second << std::endl;
	}

	// business time
	c10::Dict< c10::IValue, c10::IValue > sample_dict = sample_method( inputs, optional_inputs ).toGenericDict();

	// parse output dictionary
	torch::Tensor sample_seq_t = sample_dict.at( "S" ).toTensor();
	torch::Tensor sample_probs_t = sample_dict.at( "probs" ).toTensor();

	sample_seqs = seqs_from_tensor( sample_seq_t, pose.total_residue() );
	std::string ref_seq = pose.sequence();

	// debug outputs
	for( auto & sample_seq : sample_seqs ){
		int match = 0;
		for( core::Size ii = 0; ii < sample_seq.size(); ++ii){
			if( sample_seq[ii] == ref_seq[ii] ) match++;
		}
		TR.Debug << "Reference: " << ref_seq << std::endl;
		TR.Debug << "Sampled:   " << sample_seq << std::endl;
		TR.Debug << "Recovery = " << match / double( ref_seq.size() ) << std::endl;
	}

	return sample_seqs;
}

/// @brief Sample a single sequence using default options
std::string
ProteinMPNN::sample( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options  ){

	// the PyTorch model returns a list of strings, even when batch_size = 1
	utility::vector1< std::string > sample_seqs = sample_batch( pose, mpnn_options );

	// make sure return isn't empty
	if( sample_seqs.empty() ){
		TR.Error << "ProteinMPNN sampling failed!" << std::endl;
		// raise an exception??
		return "";
	} else {
		return sample_seqs[1];
	}
}

/// @brief Sample a single sequence using default options
std::string
ProteinMPNN::sample( const core::pose::Pose & pose ){
	TR.Debug << "Creating default ProteinMPNNOptions..." << std::endl;
	ProteinMPNNOptions mpnn_options( pose );

	return sample( pose, mpnn_options );
}

///@brief Predict amino acid probabilities/logits and fill them into maps
void
ProteinMPNN::get_probabilities_and_logits( core::pose::Pose const & pose, ProteinMPNNOptions const & mpnn_options,
                                           std::map<core::Size, utility::vector1<core::Real>>& probs_map,
                                           std::map<core::Size, utility::vector1<core::Real>>& logits_map){

    // verify options
    if (!verify_options(pose, mpnn_options)) {
        TR.Error << "Bad inputs, aborting ProteinMPNN::sample()..." << std::endl;
        throw CREATE_EXCEPTION(utility::excn::Exception,
                               "Invalid options passed to ProteinMPNN. Check log for details.");
    }

    // do we want the model's sample() or tied_sample() method
    bool tied(false);
    torch::jit::script::Method sample_method = mpnn_module.get_method("sample");
    if (!mpnn_options.tied_positions.empty()) {
        tied = true;
        sample_method = mpnn_module.get_method("tied_sample");
    }

    // prepare PyTorch inputs
    std::vector<torch::jit::IValue> inputs = generate_sample_inputs(pose, mpnn_options, tied);
    std::unordered_map<std::string, torch::jit::IValue> optional_inputs = generate_sample_optional_inputs(pose,
                                                                                                          mpnn_options,
                                                                                                          tied);
    // inference
    c10::Dict<c10::IValue, c10::IValue> sample_dict = sample_method(inputs, optional_inputs).toGenericDict();
    torch::Tensor sample_probs_t = sample_dict.at("probs").toTensor();

    // transfer to CPU if on GPU
    if (sample_probs_t.is_cuda()) {
        sample_probs_t = sample_probs_t.cpu();
    }

    // get logits from probs
    torch::Tensor sample_logits_t = torch::logit(sample_probs_t);

    // get accessor to the two tensors
    auto accessor_probs = sample_probs_t.accessor<float, 3>();
    auto accessor_logits = sample_logits_t.accessor<float, 3>();

    // iterate over tensors and fill return map
    for (int protein_pos = 0; protein_pos < accessor_probs.size(1); ++protein_pos) {
        utility::vector1<core::Real> prob_vector;
        utility::vector1<core::Real> logit_vector;

        // iterate over the third dimension (probabilities) to fill the prob_vector and logit_vector
        for (int prob_pos = 0; prob_pos < accessor_probs.size(2); ++prob_pos) {
            prob_vector.push_back(accessor_probs[0][protein_pos][prob_pos]);
            logit_vector.push_back(accessor_logits[0][protein_pos][prob_pos]);
        }
        probs_map[protein_pos + 1] = prob_vector;
        logits_map[protein_pos + 1] = logit_vector;
    }
}

/// @brief Get the citation for ProteinMPNN
/// @details TODO: fill in details for Dauparas et al.
/*static*/
basic::citation_manager::CitationCollectionBaseCOP
ProteinMPNN::get_ProteinMPNN_neural_net_citation() {
	using namespace basic::citation_manager;
	CitationCollectionOP citation(
		utility::pointer::make_shared< CitationCollection >(
		"ProteinMPNN", CitedModuleType::NeuralNetwork
		)
	);
	citation->add_citation( CitationManager::get_instance()->get_citation_by_doi( "10.1126/science.add2187" ) );
	return citation;
}


////////////////////////////////
// PRIVATE MEMBER FUNCTIONS:
////////////////////////////////

/// @brief Initialize the Pytorch model used by this class.
/// @details Called ONCE by class constructor.  Triggers read from disk!
void
ProteinMPNN::init_mpnn_model() {
	TR << "Loading TorchScript model..." << std::endl;
	try {
		mpnn_module = torch::jit::load( basic::database::full_name( "protocol_data/protein_mpnn/ProteinMPNN.pt" ) );
	} catch (const c10::Error& e) {
		TR.Error << "error loading the model\n";
		TR.Error << e.msg() << '\n';
		// raise exception??
	}
}

/// @brief Create PyTorch inputs from config options
std::vector< torch::jit::IValue >
ProteinMPNN::generate_sample_inputs( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options, bool tied ){
	std::vector<torch::jit::IValue> inputs;

	// pose-dependent inputs
	torch::Tensor coords = bb_coords_to_tensor( pose.conformation(), mpnn_options.batch_size );
	torch::Tensor randn = randn_to_tensor( pose.total_residue(), mpnn_options.batch_size );
	torch::Tensor seq = seq_to_tensor( mpnn_options.sequence, mpnn_options.batch_size );
	torch::Tensor chain_encoding = chain_encoding_to_tensor( pose, mpnn_options.batch_size );
	torch::Tensor residue_idx = residue_idx_to_tensor( pose, mpnn_options.batch_size );

	// masking inputs
	torch::Tensor coord_mask = coord_mask_to_tensor( mpnn_options.coord_mask, mpnn_options.batch_size );
	torch::Tensor chain_mask = chain_mask_to_tensor( pose, mpnn_options.chain_mask, mpnn_options.batch_size );
	//torch::Tensor chain_mask_pos = chain_mask_pos_to_tensor( pose, mpnn_options.chain_mask_pos, mpnn_options.batch_size );
	torch::Tensor pos_mask = pos_mask_to_tensor( mpnn_options.pos_mask, mpnn_options.batch_size );

	// bias inputs
	torch::Tensor omit_AAs = omit_AAs_to_tensor( mpnn_options.omit_AAs );
	torch::Tensor bias_AAs = bias_AAs_to_tensor( mpnn_options.bias_AAs );

	// tied positions
	TiedPositions tied_pos = convert_tied_positions( mpnn_options.tied_positions );

	// DO NOT CHANGE THE ORDER OF THESE INPUTS
	inputs.push_back( coords );
	inputs.push_back( randn );
	inputs.push_back( seq );
	inputs.push_back( chain_mask );
	inputs.push_back( chain_encoding );
	inputs.push_back( residue_idx );
	inputs.push_back( coord_mask );
	inputs.push_back( omit_AAs );
	inputs.push_back( bias_AAs );
	inputs.push_back( pos_mask );
	if( tied ){
		inputs.push_back( tied_pos );
	}

	return inputs;
}

/// @brief Create PyTorch optional inputs from config options
std::unordered_map< std::string, torch::jit::IValue >
ProteinMPNN::generate_sample_optional_inputs(const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options, bool tied) {
	std::unordered_map<std::string, torch::jit::IValue> optional_inputs;

	// deterministic mode
	torch::IValue deterministic( mpnn_options.deterministic_flag );
	optional_inputs["deterministic_flag"] = deterministic;

	// temperature
	float temp( mpnn_options.temperature );
	optional_inputs["temperature"] = temp;

	// per-residue omit AAs
	torch::Tensor omit_AA_mask = omit_AAs_pos_to_tensor( mpnn_options.omit_AAs_pos, mpnn_options.batch_size );
	optional_inputs["omit_AA_mask"] = omit_AA_mask;

	return optional_inputs;
}

bool
ProteinMPNN::verify_options( const core::pose::Pose & pose, const ProteinMPNNOptions & mpnn_options ){

	if( mpnn_options.batch_size != 1 ){
		TR.Error << "!!! Currently only batch_size = 1 is supported !!!" << std::endl;
		return false;
	}

	// verify sequence length and check for unsupported restypes
	if( mpnn_options.sequence.size() != pose.total_residue() ){
		TR.Error << "!!! sequence size does not match pose.total_residue() !!!" << std::endl;
		return false;
	}
	for( auto & c : mpnn_options.sequence ){
		if( AA_ALPHABET.find( c ) == std::string::npos ){
			TR.Error << "!!! Found unrecognized residue type " << c << " in input sequence !!!" << std::endl;
			return false;
		}
	}

	// verify chain_mask size
	if( mpnn_options.chain_mask.size() != pose.num_chains() ){
		TR.Error << "!!! chain_mask size does not match pose.num_chains() !!!" << std::endl;
		return false;
	}

	// verify pos_mask size
	if( mpnn_options.pos_mask.size() != pose.total_residue() ){
		TR.Error << "!!! pos_mask size does not match pose.total_residue() !!!" << std::endl;
		return false;
	}

	// verify mask size
	if( mpnn_options.coord_mask.size() != pose.total_residue() ){
		TR.Error << "!!! coord_mask size does match pose.total_residue() !!!" << std::endl;
		return false;
	}

	// verify omit_AAs
	// error if 'X' is not in omit_AAs
	// technically, this would not break the model, but I'm not sure how you'd interpret a restype prediction of 'X'...
	auto it = std::find( mpnn_options.omit_AAs.begin(), mpnn_options.omit_AAs.end(), 'X' );
	if( it == mpnn_options.omit_AAs.end() ){
		TR.Error << "!!! omit_AAs does not include restype 'X' !!!" << std::endl;
		return false;
	}

	// TODO: implement omit_AAs_pos
	// verify omit_AAs_pos
	// also verify if X included at every position?
	/*
	if( mpnn_options.omit_AAs_pos.size() != pose.total_residue() ){
		TR.Error << "!!! omit_AAs_pos size does not match pose.total_residue() !!!" << std::endl;
		return false; // TODO: implement and re-enable
	}


	for( core::Size ii = 1; ii < mpnn_options.omit_AAs_pos.size(); ++ii ){
		utility::vector1< char > omit_pos = mpnn_options.omit_AAs_pos[ ii ];

		auto it = std::find( omit_pos.begin(), omit_pos.end(), 'X' );
		if( it == omit_pos.end() ){
			TR.Error << "!!! omit_AAs_pos does not include restype 'X' at position " << ii << " !!!" << std::endl;
			return false;
		}
	}
	*/


	// verify bias_AAs size
	// remember, ProteinMPNN has 21 residue types, including unknown 'X'
	if( mpnn_options.bias_AAs.size() != AA_ALPHABET.size() ){
		TR.Error << "!!! bias_AAs size does not match the size of AA_ALPHABET (did you forget about restype 'X'?) !!!" << std::endl;
		return false;
	}

	// verify temperature
	// warning if temperature is >1.0 (completely random sampling)
	// this is allowed, but usually undesirable
	if( mpnn_options.temperature > 1.0 ){
		TR.Warning << "!!! temperature is greater than 1.0 and will sample randomly !!!" << std::endl;
	}

	// verify tied_pos
	// just make sure residue indices exist
	for( auto & tied_set : mpnn_options.tied_positions ){
		for( auto & res : tied_set ){
			if( res < 1 || res > pose.total_residue() ){
				TR.Error << "!!! tied_positions includes invalid res index " << res << " for pose of size " << pose.total_residue() << " !!!" << std::endl;
			}
		}
	}

	return true;
}

ProteinMPNNOptions::ProteinMPNNOptions( const core::pose::Pose & pose ){
	reset( pose );
}

void
ProteinMPNNOptions::reset( const core::pose::Pose & pose ){

	// default sequence: pose sequence
	// replace unrecognized residues with 'X'
	sequence = "";
	for( auto & c : pose.sequence() ){
		if( AA_ALPHABET.find( c ) == std::string::npos ){
			sequence += 'X';
		} else {
			sequence += c;
		}
	}

	// default mask: include all residue coords
	coord_mask.clear();
	coord_mask.resize( pose.total_residue(), true );

	// default chain_mask: predict everything
	// NB: chain_mask lists chain_ID (NOT chain lettering)
	chain_mask.clear();
	chain_mask.resize( pose.num_chains(), true );

	// default pos_mask: predict all residues
	pos_mask.clear();
	pos_mask.resize( pose.total_residue(), true );

	// default omit_AAs: omit restype 'X'
	omit_AAs.clear();
	omit_AAs.push_back( 'X' );

	// default omit_AAs_pos: all residues omit unknown restype 'X'
	omit_AAs_pos.clear();
	for( core::Size ii = 1; ii <= pose.total_residue(); ++ii ){
		utility::vector1< char > res_omit_AAs;
		res_omit_AAs.push_back( 'X' );
		omit_AAs_pos.push_back( res_omit_AAs );
	}

	// default bias_AAs: zeroes
	bias_AAs.clear();
	bias_AAs.resize( AA_ALPHABET.size(), 0.0 );

	// default temperature: 0.1
	temperature = 0.1;

	// default tied_positions: none tied, all positions independent
	tied_positions.clear();

	// default batch_size: 1
	batch_size = 1;

	// ... mask coords and positions for non-protein residues
	for( core::Size ii = 1; ii <= pose.total_residue(); ++ii ){
		if( !pose.residue_type( ii ).is_protein() ){
			coord_mask[ ii ] = false;
			pos_mask[ ii ] = false;
		}
	}

	// default deterministic_flag: false
	// NB: setting this to true will override options like temperature
	deterministic_flag = false;
}

} //protein_mpnn
} //protocols

#endif //USE_TORCH

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.cc
/// @brief A class for predicting amino acid probabilities P(Sequence|Sequence) using the ESM language model.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adapted from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)

// project headers
#include <protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.hh>

// // utility headers
#include <utility>
#include <utility/pointer/owning_ptr.hh>
#include <utility/file/file_sys_util.hh>

// basic headers
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/execute.hh>
#include <basic/tensorflow_manager/RosettaTensorflowManager.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/machine_learning.OptionKeys.gen.hh>

// core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>

// protocol headers

// STL headers:
#include <chrono>
#include <iostream>
#include <cmath>
#include <boost/format.hpp>

static basic::Tracer TR( "protocols.esm_perplexity.EsmPerplexityTensorflowProtocol" );


namespace protocols {
namespace esm_perplexity {

// defines the base directory where the models are stored on the MeilerLab IWE GitLab
std::string const EsmPerplexityTensorflowProtocol::base_url_ = "https://git.iwe-lab.de/moritzertelt/ML_graphs/-/archive/main/ML_graphs-main.tar.gz?path=tensorflow_graphs/ESM/";
// defines which models are available
std::vector< std::string > const EsmPerplexityTensorflowProtocol::available_models_ = {"esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D", "esm2_t30_150M_UR50D", "esm2_t33_650M_UR50D"};
// define the ESM alphabet used for tokenization, first 4 letters are placeholders for special characters
std::string const EsmPerplexityTensorflowProtocol::alphabet_ = "XXXXLAGVSERTIDPKQNFYMHWC";

/// @brief Constructor using the model name to initialize the particular ESM model
EsmPerplexityTensorflowProtocol::EsmPerplexityTensorflowProtocol(std::string const & model) : model_name_( model ) {
	set_defaults();
#ifdef USE_TENSORFLOW
    set_tensorflow_session(get_tensorflow_session());
#endif //USE_TENSORFLOW
}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
basic::tensorflow_manager::RosettaTensorflowProtocolBaseOP
EsmPerplexityTensorflowProtocol::clone() const {
	return utility::pointer::make_shared<EsmPerplexityTensorflowProtocol>(*this);
}

/// @brief Get the name of this protocol.
std::string
EsmPerplexityTensorflowProtocol::name() const {
	return "EsmPerplexityTensorflowProtocol";
}

/// @brief Checks if the given model name is one of the available models
bool
EsmPerplexityTensorflowProtocol::is_model_valid( std::string const & model_name) {
	return std::find(available_models_.begin(), available_models_.end(), model_name) != available_models_.end();
}

/// @brief Returns a string containing the names of all available models, separated by commas
std::string
EsmPerplexityTensorflowProtocol::concatenate_model_names() {
	std::string model_names;
	for ( size_t i = 0; i < available_models_.size(); ++i ) {
		model_names += available_models_[i];
		if ( i != available_models_.size() - 1 ) {
			model_names += ", ";
		}
	}
	return model_names;
}

void
EsmPerplexityTensorflowProtocol::set_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	auto_download( basic::options::option[ OptionKeys::machine_learning::auto_download ]());
}

/// @brief Set the auto_download bool
void
EsmPerplexityTensorflowProtocol::auto_download( bool setting ) { auto_download_ = setting; }

/// @brief Given a selection, return the number of selected residues.
utility::vector1<core::Size>
EsmPerplexityTensorflowProtocol::get_selected_residues(
	core::select::residue_selector::ResidueSubset const &selected
) {
	utility::vector1<core::Size> returnvec;
	for ( core::Size i(1), imax(selected.size()); i <= imax; ++i ) {
		if ( selected[i] ) returnvec.push_back(i);
	}
	return returnvec;
}

#ifdef USE_TENSORFLOW
/// @brief Set the tensorflow session
void EsmPerplexityTensorflowProtocol::set_tensorflow_session(
        basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session_in) {
    debug_assert(session_in != nullptr);
    tensorflow_session_ = std::move(session_in);
}

/// @brief Get the Tensorflow session used by this class.
/// @details Implemented by derived classes.
basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP
EsmPerplexityTensorflowProtocol::get_tensorflow_session(
) const {
    // check if the provided model name matches one of the available ESM models
    if (!is_model_valid(model_name_)) {
        utility_exit_with_message("Invalid model name: " + model_name_ + ". Please provide one of the following models: " + concatenate_model_names());
    }
    std::string const path_to_model = "protocol_data/tensorflow_graphs/tensorflow_graph_repo_submodule/ESM/" + model_name_;
    // if the flag is set attempt to automatically download the models if they are missing
    // else print out instructions to either set the flag or download everything manually
    // check if the model files are existing, and if not attempt to download them
    download_model_if_not_existing( path_to_model, base_url_, auto_download_ );

    return basic::tensorflow_manager::RosettaTensorflowManager::get_instance()->get_session(path_to_model,
                                                                                            "serve");
}

/// @brief Analyze a given position in a pose, returning the probability of all possible amino acids.
std::map<core::Size, utility::vector1<core::Real>>
EsmPerplexityTensorflowProtocol::compute_logits(
        core::pose::Pose const &pose,
        core::select::residue_selector::ResidueSubset const &selection,
        core::select::residue_selector::ResidueSubset  const &attention_mask_selection,
        bool const &multirun
) const {

    // Selected pose residues:
    utility::vector1<core::Size> const selected_residues(get_selected_residues(selection));
    // Selected residues to be masked
    utility::vector1<core::Size> const attention_mask_vec(get_selected_residues(attention_mask_selection));

    // get pose sequence and length
    std::string const pose_sequence = pose.sequence();

    // setting up input and output to run all predictions at once, instead of one at a time
    utility::vector1<basic::tensorflow_manager::RosettaTensorflowTensorContainer<int> > input_tensor_mask_vec(
            selected_residues.size());
    utility::vector1<basic::tensorflow_manager::RosettaTensorflowTensorContainer<int> > input_tensor_seq_vec(
            selected_residues.size());
    utility::vector1<basic::tensorflow_manager::RosettaTensorflowTensorContainer<float> > output_tensor_vec(
            selected_residues.size());
//	for ( auto const position : selected_residues ) { // position is 1 indexed
    for (size_t position_index = 1; position_index <= selected_residues.size(); ++position_index) {

        core::Size position = selected_residues[position_index];
        // check that the selected position is a canonical amino acid,
        // we can handle NCAAs being present in the pose but they should not be part of the selection
        runtime_assert_string_msg(pose.residue(position).has_property("CANONICAL_AA"),
                                  "Can only calculate logits for residues that are canonical amino acids!");
        // get chain of selected residue
        core::Size const chain = pose.residue(position).chain();
        // get chain sequence
        std::string const chain_sequence = pose.chain_sequence(chain);
        // get chain length
        long const chain_length = static_cast<long>(chain_sequence.size());

        // Creating a tensorflow container for the input & output tensor
        // Size is the length of the sequence + 2 (Start and End token)
        input_tensor_seq_vec[position_index].initialize(TF_INT32, {1, chain_length + 2});
        input_tensor_mask_vec[position_index].initialize(TF_INT32, {1, chain_length + 2});
        output_tensor_vec[position_index].initialize(TF_FLOAT, {1, chain_length + 2, 33}, 0.0);

        copy_feat_to_tensor(pose, position, attention_mask_vec, input_tensor_seq_vec[position_index],
                            input_tensor_mask_vec[position_index]);
    }
    // concatenating the tensor vectors for the two inputs
    utility::vector1<utility::vector1<basic::tensorflow_manager::RosettaTensorflowTensorContainer<int> > > multirun_input_vec;
    multirun_input_vec.push_back(input_tensor_mask_vec);
    multirun_input_vec.push_back(input_tensor_seq_vec);
    TR << "================== Starting prediction ===============" << std::endl;
    // setup the input names of the tensorflow graph
    utility::vector1<std::string> input_names;
    input_names.push_back("serving_default_attention_mask");
    input_names.push_back("serving_default_input_ids");
    // setup name for output tensor of TF graph
    std::string const output_name = "StatefulPartitionedCall";
    // setup session
    basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session = EsmPerplexityTensorflowProtocol::tensorflow_session();

    if (multirun) {
        // Calculate logits for all selected positions simultaneously
        // Chrono timing required by run_session
        std::chrono::duration<double, std::micro> runtime{};

        // predicting based on the input tensor
        session->multirun_multiinput_session(
                input_names,
                output_name,
                multirun_input_vec,
                output_tensor_vec,
                runtime);

#ifdef USE_TENSORFLOW_GPU
        TR << "Computed logits using Tensorflow-GPU in " << runtime.count() << " microseconds." << std::endl;
#else //USE_TENSORFLOW_GPU
#ifdef USE_TENSORFLOW_CPU
        TR << "Computed logits using Tensorflow-CPU in " << runtime.count() << " microseconds." << std::endl;
#else
        utility_exit_with_message( "Either USE_TENSORFLOW_GPU or USE_TENSORFLOW_CPU should be defined! If not you are probably in deep trouble." );
#endif //USE_TENSORFLOW_CPU
#endif //USE_TENSORFLOW_GPU

    } else {

        // calculate logits for each position one by one
        for (size_t position_index = 1; position_index <= selected_residues.size(); ++position_index) {

            // getting the right tensors from all input tensors for each position in selected_residues
            utility::vector1<basic::tensorflow_manager::RosettaTensorflowTensorContainer<int> > singlerun_input_vec;
            singlerun_input_vec.push_back(multirun_input_vec[1][position_index]);
            singlerun_input_vec.push_back(multirun_input_vec[2][position_index]);

            // Chrono timing required by run_session
            std::chrono::duration<double, std::micro> runtime{};

            session->run_multiinput_session(
                    input_names,
                    output_name,
                    singlerun_input_vec,
                    output_tensor_vec[position_index], // use relevant output tensor from vector
                    runtime);

#ifdef USE_TENSORFLOW_GPU
            TR << "Computed logits using Tensorflow-GPU in " << runtime.count() << " microseconds." << std::endl;
#else //USE_TENSORFLOW_GPU
#ifdef USE_TENSORFLOW_CPU
            TR << "Computed logits using Tensorflow-CPU in " << runtime.count() << " microseconds." << std::endl;
#else
            utility_exit_with_message( "Either USE_TENSORFLOW_GPU or USE_TENSORFLOW_CPU should be defined! If not you are probably in deep trouble." );
#endif //USE_TENSORFLOW_CPU
#endif //USE_TENSORFLOW_GPU
        }
    }
    // go through positions
    std::map<core::Size, utility::vector1<core::Real>> logit_map;
    for (size_t position_index = 1; position_index <= selected_residues.size(); ++position_index) {

        utility::vector1<core::Real> logit_vec;
        // get position
        core::Size position = selected_residues[position_index];
        // only get the logits from the masked position
        for (size_t logit_index = 1; logit_index < 34; ++logit_index) {
            // fill the vector for calculating log_softmax later
            core::Real const logit = output_tensor_vec[position_index](1, position + 1, logit_index);
            logit_vec.push_back(logit);
        }
        logit_map[position] = logit_vec;
    }
    return logit_map;
}

#endif //USE_TENSORFLOW

#ifdef USE_TENSORFLOW

/// @brief Calculate features needed for ESM prediction and copy them to input tensors.
void
EsmPerplexityTensorflowProtocol::copy_feat_to_tensor(
        core::pose::Pose const &pose,
        core::Size const &position,
        utility::vector1<core::Size> const &attention_mask_vec,
        basic::tensorflow_manager::RosettaTensorflowTensorContainer<int> &input_tensor_seq,
        basic::tensorflow_manager::RosettaTensorflowTensorContainer<int> &input_tensor_mask
) {
    // get chain of selected residue
    core::Size const chain = pose.residue(position).chain();
    // get chain sequence
    std::string const chain_sequence = pose.chain_sequence(chain);
    // get start of chain
    core::Size const chain_begin = pose.chain_begin(chain);

    // define an AA alphabet, based on the ESM model
    // <mask> token has index 32, start (<cls>) is 0 and end (<eos>) is 2
    std::vector<std::string> alphabet_vector = {"<cls", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E",
                                                "R", "T", "I", "D", "P", "K", "Q", "N", "F", "Y", "M", "H", "W",
                                                "C", "X", "B", "U", "Z", "O", ".", "-", "<null_1>", "<mask>"};
    core::Size const chain_length = chain_sequence.size();
    // Reminder: rosetta tensorflow containers are 1 indexed
    // set the first value of the tensor to the start token
    input_tensor_seq(1) = 0;

    // set the last value of the tensor to the end of line token
    input_tensor_seq(chain_length + 2) = 2;
    // fill the rest depending on the sequence
    for (core::Size i = 0; i < chain_length; ++i) {
        // get amino acid identity and input into the tensor
        // check if its a canonical AA to avoid mixups, else set to ? which
        // will be converted to the <unknown> token later
        char aa_name;
        if (pose.residue(i + chain_begin).has_property("CANONICAL_AA")) {
            aa_name = chain_sequence[i];
        } else {
            aa_name = '?';
            TR.Warning << "Residue at position " << i + chain_begin
               << " is not a canonical AA and token will be set to unknown but you might want to mask it instead." << std::endl;
        }
        // mask the position that is being predicted
        if (i == ( position - chain_begin )) {
            input_tensor_seq(i + 2) = 32;
        }
            // tokenize the rest of the sequence
        else {
            // need an int here since the Tensorflow tensor expects one
            int const esm_alphabet_index =
                    static_cast<int>(std::find(alphabet_vector.begin(), alphabet_vector.end(), std::string(1, aa_name)) -
                    alphabet_vector.begin());
            if (esm_alphabet_index == 33) {
                input_tensor_seq(i + 2) = 3; // if amino acid isn't in alphabet use unknown token
            } else {
                input_tensor_seq(i + 2) = esm_alphabet_index;
            }
        }
        // setup the masking tensor, used if we'd want to ignore parts of the sequence
        // here I just make sure its properly setup, it gets modified with the residue_selector later
        input_tensor_mask(i + 2) = 1;
    }
    // don't forget first and last value of mask tensor
    input_tensor_mask(1) = 1;
    input_tensor_mask(chain_length + 2) = 1;

    // loop through the attention mask selection provided by the user
    // check whether the selected residue is on the same chain as the position to be predicted and mask
    for (auto const &mask_position : attention_mask_vec ) {
        core::Size const mask_chain = pose.residue( mask_position ).chain();
        if ( mask_chain == chain ) {
            input_tensor_mask( mask_position - chain_begin + 1 ) = 0;
        }
    }
}
#endif // USE_TENSORFLOW

/// @brief Function to return the softmax of logits resulting in probabilities
void
EsmPerplexityTensorflowProtocol::softmax(
	std::map< core::Size, utility::vector1< core::Real >> const &logit_map,
	std::map< core::Size, utility::vector1< core::Real >> &softmax_map
) {
	for ( const auto& pair : logit_map ) {
		core::Size selected_residue = pair.first;
		utility::vector1< core::Real > logit_vec = pair.second;

		core::Real max_val = *std::max_element(logit_vec.begin(), logit_vec.end());

		for ( auto& logit : logit_vec ) {
			logit -= max_val;
		}

		core::Real sum_exp = 0.0;
		for ( const auto& scaled_logit : logit_vec ) {
			sum_exp += std::exp(scaled_logit);
		}

		utility::vector1< core::Real > softmax_vec;
		for ( const auto& scaled_logit : logit_vec ) {
			softmax_vec.push_back(std::exp(scaled_logit) / sum_exp);
		}

		softmax_map[selected_residue] = softmax_vec;
	}

}

/// @brief Downloads model from GitLab if the specified path does not exist or is missing crucial files
void
EsmPerplexityTensorflowProtocol::download_model_if_not_existing( std::string const & path_to_model, std::string const & base_url, bool const & auto_download ) {
	// Form GitLab URL
	std::size_t model_name_start = path_to_model.rfind('/') + 1;
	std::string model_name = path_to_model.substr(model_name_start);
	std::string gitlab_url = base_url + model_name;
	// get the full path to the model necessary for later commands
	std::string const full_path_to_model = basic::database::full_name( path_to_model );
	// check if all the needed files are present, this is specific to tensorflow
	utility::vector1<std::string> files;
	utility::file::list_dir(full_path_to_model, files);
	// file_missing bool is true if any of the following files are missing
	bool files_are_missing = (std::find(files.begin(), files.end(), "variables") == files.end()) ||
		(std::find(files.begin(), files.end(), "keras_metadata.pb") == files.end()) ||
		(std::find(files.begin(), files.end(), "saved_model.pb") == files.end());
	// Check for directory existence and has all necessary files
	if ( !utility::file::file_exists(full_path_to_model ) || files_are_missing  ) {
		// automatically download if the user has set the -auto_download flag
		if ( auto_download ) {
			std::string const tar_file = model_name + ".tar.gz";
			TR.Info << "Downloading missing model files.... " << std::endl;
			// TODO: remove certificate flag after updating ssl stuff on our gitlab
			std::string message = "Downloading model from GitLab";
			std::string command = "wget";
			std::vector<std::string> args = {gitlab_url, "--progress=bar:force", "--no-check-certificate", "-O", tar_file};
			basic::ExecutionResult result = basic::execute(message, command, args, false, true );
			if ( result.result != 0 ) {
				utility_exit_with_message("Download failed. Please manually download the model from: " + gitlab_url + " , than extract the model directory to " + full_path_to_model );
			}
			TR.Info << "Model successfully downloaded." << std::endl;
			message = "Extracting model files";
			command = "tar";
			args = {"-xf", tar_file};
			result = basic::execute(message, command, args, false, false );
			if ( result.result != 0 ) {
				utility_exit_with_message( "Extraction failed. Please manually extract the model from the downloaded tar.gz file using: tar -xf " + tar_file );
			}
			// Move the files to the correct location
			std::string src_dir = "ML_graphs-main-tensorflow_graphs-ESM-" + model_name + "/tensorflow_graphs/ESM/" + model_name ;
			message = "Moving model files";
			command = "mv";
			args = {src_dir, full_path_to_model};
			result = basic::execute(message, command, args, false, false );
			if ( result.result != 0 ) {
				utility_exit_with_message( "Could not move the model to the appropriate folder. Please move or copy " + src_dir + " to " + full_path_to_model);
			}
			// Remove the leftover directories
			std::string remove_dir = "ML_graphs-main-tensorflow_graphs-ESM-" + model_name;
			message = "Removing unnecessary directories";
			command = "rm";
			args = {"-r", remove_dir};
			result = basic::execute(message, command, args, false, true);
			if ( result.result != 0 ) {
				TR.Warning << "Could not remove unnecessary folders created from model extraction, you might want to delete " + remove_dir << std::endl;
			}
		} else {
			utility_exit_with_message("The required model is missing. You can either automatically download and extract it by setting the -auto_download flag or manually download the model from: " + gitlab_url + " and than extract it to " + full_path_to_model );
		}
	}
}

/// @brief Get the citation for ESM
/*static*/
basic::citation_manager::CitationCollectionBaseCOP
EsmPerplexityTensorflowProtocol::get_ESM_neural_net_citation() {
	using namespace basic::citation_manager;
	CitationCollectionOP citation(
		utility::pointer::make_shared< CitationCollection >(
		"ESM", CitedModuleType::NeuralNetwork
		)
	);
	citation->add_citation( CitationManager::get_instance()->get_citation_by_doi( "10.1126/science.ade2574" ) );
	return citation;
}


} // esm_perplexity
} // protocols

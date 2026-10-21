// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/inverse_folding/MIFST.cc
/// @brief A class for using the MIF-ST model developed by Yang et al.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#ifdef USE_TORCH

// Project headers:
#include <protocols/inverse_folding/MIFST.hh>

// Basic headers:
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/machine_learning.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <basic/execute.hh>

// core headers
#include <core/conformation/Residue.hh>
#include <core/chemical/Atom.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>

// numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

// Utility headers:
#include <utility/pointer/memory.hh>
#include <utility/file/file_sys_util.hh>

// torch headers:
#include <torch/csrc/autograd/generated/variable_factories.h>
#include <torch/script.h>

// protocol headers:
#include <protocols/inverse_folding/util.hh>


static basic::Tracer TR( "protocols.inverse_folding.MIFST" );


namespace protocols {
namespace inverse_folding {

std::string const MIFST::GITLAB_URL_ = "https://git.iwe-lab.de/moritzertelt/ML_graphs/-/archive/main/ML_graphs-main.tar.gz?path=pytorch_graphs/MIF-ST";

/// @brief constructor.
MIFST::MIFST() {
    set_defaults();
    init_mifst_model();
}

////////////////////////////////
// PUBLIC MEMBER FUNCTIONS:
////////////////////////////////
std::map< core::Size, std::map< core::chemical::AA, core::Real > >
MIFST::sample( core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & residue_selection, core::select::residue_selector::ResidueSubset const & feature_selection, bool const multirun, bool const use_gpu ){

    // get the logits (runs inference)
    torch::Tensor logits = predict( pose, residue_selection, feature_selection, multirun, use_gpu );

    // Apply softmax to transform logits into probabilities
    torch::Tensor probabilities = torch::softmax(logits, -1);

    // Populate the predictions map
    std::map<core::Size, std::map<core::chemical::AA, core::Real>> predictions;
    int64_t torch_index = 0;
    int64_t feature_index = 0;
    for (core::Size residue_pos = 1; residue_pos <= residue_selection.size(); ++residue_pos) {
        if (feature_selection[residue_pos]) {
            if (residue_selection[residue_pos]) {
                std::map<core::chemical::AA, core::Real> aa_prob_map;
                for (core::Size aa_pos = 0; aa_pos < 20; ++aa_pos) { // first 20 positions are the canonicals
                    core::chemical::AA aa = core::chemical::aa_from_oneletter_code(protocols::inverse_folding::PROTEIN_ALPHABET[aa_pos]);
                    auto prob = probabilities[torch_index][feature_index][static_cast<int64_t>(aa_pos)].item<core::Real>();
                    aa_prob_map[aa] = prob;
                }
                predictions[residue_pos] = aa_prob_map;
                torch_index += 1;
            }
            feature_index += 1;
        }
    }
    return predictions;
}

torch::Tensor
MIFST::predict( core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & residue_selection, core::select::residue_selector::ResidueSubset const & feature_selection, bool const multirun, bool const use_gpu ){

    // checking whether inference can and should be run on GPU
    bool cuda_wanted_and_available = false;
    if ( use_gpu ) {
#ifndef USE_TORCH_GPU
       TR.Warning << "If you'd like to run PyTorch on the GPU you need to download the appropriate libtorch library and compile with extras=pytorch_gpu! Falling back to CPU instead." << std::endl;
#else
       if ( torch::hasCUDA() ) {
          cuda_wanted_and_available = true;
       } else {
           TR.Warning << "You need to download the appropriate libtorch library to run inference on the GPU! Falling back to CPU instead." << std::endl;
       }
#endif //USE_TORCH_GPU
    }

    // Process coordinates and get features
    auto result = protocols::inverse_folding::process_coords(pose, feature_selection);
    torch::Tensor dist, omega, theta, phi;
    std::tie(dist, omega, theta, phi) = result;

    // setup the sequence and mask vector
    std::vector<std::string> sequences;
    std::vector<core::Size> mask_positions;

    //  creating the sequence based on feature_selection, while checking for NCAAs
    std::string feature_sequence;
    for (core::Size i = 1; i <= pose.total_residue(); ++i) {
        if (feature_selection[i]) {
            core::chemical::ResidueType const & res_type = pose.residue_type(i);
            if (res_type.is_canonical()) {
                feature_sequence += res_type.name1();
            } else {
                feature_sequence += protocols::inverse_folding::MASK; // if its not a canonical we mask that position as well
            }
        }
    }

    // Create mask_positions based on residue_selection
    for (core::Size residue_index = 1; residue_index <= residue_selection.size(); ++residue_index) {
        if (residue_selection[residue_index]) {
            mask_positions.push_back(residue_index - 1); // -1 here since I don't want to later switch when handling torch tensors
            sequences.push_back(feature_sequence); // one sequence for each masked_position
        }
    }

    // Collate data
    auto collated_data = protocols::inverse_folding::collate_structure(sequences,
                                                             dist,
                                                             omega,
                                                             theta,
                                                             phi,
                                                             mask_positions);

    torch::Tensor src, nodes, edges, connections, edge_mask;
    std::tie(src, nodes, edges, connections, edge_mask) = collated_data;

    torch::Tensor logits;
    // Move inputs and model to GPU if CUDA is available and wanted
    if (cuda_wanted_and_available) {
        mifst_module_.to(at::Device(at::kCUDA));
	src = src.to(torch::device(torch::kCUDA));
	nodes = nodes.to(torch::device(torch::kCUDA));
	edges = edges.to(torch::device(torch::kCUDA));
	connections = connections.to(torch::device(torch::kCUDA));
	edge_mask = edge_mask.to(torch::device(torch::kCUDA));
    }

    // decide whether we run inference on all residues of the selection at once, or one-by-one
    {
	    // enabling no_grad 	    
	    torch::NoGradGuard no_grad;
	    TR << "Starting prediction..." << std::endl;
	    if (multirun) {

		// Prepare inputs for the model
		std::vector<torch::jit::IValue> inputs;
		inputs.emplace_back(src);
		inputs.emplace_back(nodes);
		inputs.emplace_back(edges);
		inputs.emplace_back(connections);
		inputs.emplace_back(edge_mask);
		// inference time
		logits = mifst_module_(inputs).toTensor();

		// Move logits back to CPU if they were on GPU
		if (cuda_wanted_and_available) {
		    logits = logits.to(torch::device(torch::kCPU));
		}
		// Clear IValues explicitly to free up GPU memory
		inputs.clear();


	    } else {
		std::vector<torch::Tensor> logits_list;
		for (int64_t i = 0; i < src.size(0); ++i) { //int64_t because we are in torch land
		    // Prepare inputs for the model for each sequence
		    std::vector<torch::jit::IValue> inputs;
		    inputs.emplace_back(src[i].unsqueeze(0));
		    inputs.emplace_back(nodes[i].unsqueeze(0));
		    inputs.emplace_back(edges[i].unsqueeze(0));
		    inputs.emplace_back(connections[i].unsqueeze(0));
		    inputs.emplace_back(edge_mask[i].unsqueeze(0));

		    // inference time for each sequence
		    torch::Tensor single_logit = mifst_module_(inputs).toTensor();

		    // Move logits back to CPU if they were on GPU
		    if (cuda_wanted_and_available) {
			single_logit = single_logit.to(torch::device(torch::kCPU));
		    }

		    logits_list.push_back(single_logit.squeeze(0));
		    // Clear IValues explicitly to free up GPU memory
		    inputs.clear();
		}
		// Stack the logits to form the final tensor
		logits = torch::stack(logits_list, 0);
	    }
	    TR << "... finished!" << std::endl;
    }

    return logits;
}

void
MIFST::set_defaults() {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    auto_download( basic::options::option[ OptionKeys::machine_learning::auto_download ]());
}

/// @brief Set the auto_download bool
void
MIFST::auto_download( bool setting ) { auto_download_ = setting; }

/// @brief Get the citation for MIF-ST
/*static*/
basic::citation_manager::CitationCollectionBaseCOP
MIFST::get_MIFST_neural_net_citation() {
    using namespace basic::citation_manager;
    CitationCollectionOP citation(
            utility::pointer::make_shared< CitationCollection >(
                    "MIFST", CitedModuleType::NeuralNetwork
            )
    );
    citation->add_citation( CitationManager::get_instance()->get_citation_by_doi( "10.1101/2022.05.25.493516" ) );
    return citation;
}

////////////////////////////////
// PRIVATE MEMBER FUNCTIONS:
////////////////////////////////

/// @brief Initialize the Pytorch model used by this class.
/// @details Called ONCE by class constructor.  Triggers read from disk!
void
MIFST::init_mifst_model() {
    TR << "Loading TorchScript model..." << std::endl;
    // if the auto_download flag is set attempt to automatically download the models if they are missing
    // else print out instructions to either set the flag or download everything manually
    // check if the model files are existing, and if not attempt to download them
    std::string const path_to_model = "protocol_data/inverse_folding/mifst/";
    download_model_if_not_existing( path_to_model, auto_download_ );
    try {
        mifst_module_ = torch::jit::load( basic::database::full_name( "protocol_data/inverse_folding/mifst/MIF-ST_traced.pt" ) );
        mifst_module_.eval();
        TR << "... Success!" << std::endl;
    } catch (const c10::Error& e) {
        TR.Error << "error loading the model\n";
        TR.Error << e.msg() << '\n';
    }
}

/// @brief Downloads model from GitLab if the specified path does not exist or is missing crucial files
void
MIFST::download_model_if_not_existing( std::string const & path_to_model, bool const & auto_download ) {
    // get the full path to the model necessary for later commands
    std::string const full_path_to_model = basic::database::full_name( path_to_model );
    // check if all the needed files are present, this is specific to tensorflow
    utility::vector1<std::string> files;
    utility::file::list_dir(full_path_to_model, files);
    // file_missing bool is true if any of the following files are missing
    bool files_are_missing = (std::find(files.begin(), files.end(), "MIF-ST_traced.pt") == files.end()) ||
            (std::find(files.begin(), files.end(), "LICENSE") == files.end());

    // Check for directory existence and has all necessary files
    if ( !utility::file::file_exists(full_path_to_model ) || files_are_missing  ) {
        // automatically download if the user has set the -auto_download flag
        if ( auto_download ) {
            std::string const tar_file = "MIF-ST.tar.gz";
            TR.Info << "Downloading missing model files.... " << std::endl;
            // TODO: remove certificate flag after updating ssl stuff on our gitlab
            std::string message = "Downloading model from GitLab, only happens once at first use (2.5GB in total).";
            std::string command = "wget";
            std::vector<std::string> args = {GITLAB_URL_, "--progress=bar:force", "--no-check-certificate", "-O", tar_file};
            basic::ExecutionResult result = basic::execute(message, command, args, false, true );
            if ( result.result != 0 ) {
                utility_exit_with_message("Download failed. Please manually download the model from: " + GITLAB_URL_ + " , than extract the model directory to " + full_path_to_model );
            }
            TR.Info << "Model successfully downloaded." << std::endl;
            message = "Extracting model files";
            command = "tar";
            args = {"-xf", tar_file, "--strip-components=3", "-C", full_path_to_model };
            result = basic::execute(message, command, args, false, false );
            if ( result.result != 0 ) {
                utility_exit_with_message( "Extraction failed. Please manually extract the model from the downloaded tar.gz file using: tar -xf " + tar_file + " and move it to " + full_path_to_model);
            }
        } else {
            utility_exit_with_message("The required model is missing. You can either automatically download and extract it by setting the -auto_download flag or manually download the model from: " + GITLAB_URL_ + " and than extract it to " + full_path_to_model );
        }
    }
}


} //inverse_folding
} //protocols

#endif //USE_TORCH

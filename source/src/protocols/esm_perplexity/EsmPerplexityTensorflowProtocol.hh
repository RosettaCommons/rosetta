// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.hh
/// @brief A class for predicting probabilities of amino acids P(Sequence|Sequence) using the ESM language model.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adapted from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)
#ifndef INCLUDED_protocols_esm_perplexity_EsmPerplexityTensorflowProtocol_hh
#define INCLUDED_protocols_esm_perplexity_EsmPerplexityTensorflowProtocol_hh


#ifdef USE_TENSORFLOW
#include <tensorflow/c/c_api.h>
#endif

// protocol headers
#include <protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/io/ozstream.hh>

// core headers
#include <core/sequence/SequenceProfile.hh>

// Basic headers
#include <basic/citation_manager/CitationCollectionBase.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowProtocolBase.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh>

// Core headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// C++ header
#include <map>

#ifdef    SERIALIZATION

// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace esm_perplexity {

/// @brief Protocol to predict amino acid probabilities using the ESM language model family
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)
class EsmPerplexityTensorflowProtocol : public basic::tensorflow_manager::RosettaTensorflowProtocolBase {

public:

	/// @brief Default constructor.
	EsmPerplexityTensorflowProtocol() = default;

	/// @brief Copy constructor.
	EsmPerplexityTensorflowProtocol(EsmPerplexityTensorflowProtocol const &) = default;

	/// @brief Destructor
	~EsmPerplexityTensorflowProtocol() override= default;

	/// @brief Get the name of this protocol.
	std::string name() const override;

	/// @brief path to the model
	std::string model_name_;

	/// @brief The ESM alphabet used for tokenization
	static std::string const alphabet_;

	/// @brief constructor with name of the model
	explicit EsmPerplexityTensorflowProtocol( std::string const & model );

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	basic::tensorflow_manager::RosettaTensorflowProtocolBaseOP clone() const override;

	// @brief get the citation for ESM
	static
	basic::citation_manager::CitationCollectionBaseCOP
	get_ESM_neural_net_citation();

	/// @brief The tensorflow session
#ifdef USE_TENSORFLOW
	basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP tensorflow_session_;
#endif // USE_TENSORFLOW

public:
#ifdef USE_TENSORFLOW
	/**
	* @brief Calculate features needed for ESM prediction and copy them to input tensors.
	* @details Given a pose and an already-allocated (but empty) input tensor, store the relevant pose data for the
	ESM prediction in the tensor.
	* @param[in] pose
	* @param[in] position The position of the residue to predict logits for.
	* @param[in] attention_mask_vec A vector specifying which residues will be masked (aka ignored by the model prediction).
	* @param[in, out] input_tensor_seq The input tensor which will be filled with the tokenized sequence.
	* @param[in, out] input_tensor_mask The second input tensor which will be modified according to the attention_mask_vec.
	* @remark This method is public to make it accessible to unit tests.
	*/
	static void
	copy_feat_to_tensor(
		core::pose::Pose const & pose,
		core::Size const & position,
		utility::vector1<core::Size> const & attention_mask_vec,
		basic::tensorflow_manager::RosettaTensorflowTensorContainer< int > & input_tensor_seq,
		basic::tensorflow_manager::RosettaTensorflowTensorContainer< int > & input_tensor_mask
	);
#endif // USE_TENSORFLOW

	/// @brief Apply softmax function to logits
	/// @param[in] logit_map A map of logits, where each key represents a residue and its predicted logits
	/// @param[in, out] softmax_map A map of softmax, where each key represents a residue and its predicted probabilities.
	static void softmax(std::map< core::Size, utility::vector1< core::Real >> const &logit_map,
		std::map< core::Size, utility::vector1< core::Real >> &softmax_map);
	/// @brief set default options
	void set_defaults();
	/// @brief get auto_download value
	bool auto_download() const { return auto_download_; }
	/// @brief set the auto_download value
	void auto_download( bool setting );


#ifdef USE_TENSORFLOW
    /// @brief Predict logits using ESM.
    /// @details Uses the ESM language model to predict amino acid probabilities (logits) given P(Sequence|Sequence) for a selection of residues.
    /// @param[in] pose The pose attributed to the residues to predict on.
    /// @param[in] selection The residue selection to predict logits for.
    /// @param[in] attention_mask_selection A residue selector specifying which residues get masked (basically ignored by the model).
    /// @param[in] multirun Bool that can be set to run predictions for all positions simultaneously.
    /// @returns A map where each key represents a residue from the residue selector mapped to a vector of logits
    /// @remark The returned map still contains the logits of the special characters (e.g. "cls" or "eos" token)
    std::map< core::Size, utility::vector1< core::Real >> compute_logits(
            core::pose::Pose const &pose,
            core::select::residue_selector::ResidueSubset const & selection,
            core::select::residue_selector::ResidueSubset const & attention_mask_selection,
            bool const & multirun	) const ;

     /// @brief Get the Tensorflow session used by this class.
    basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP get_tensorflow_session( ) const ;

    /// @brief Allow access to the Tensorflow session
    inline basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP tensorflow_session() const { return tensorflow_session_;}
#endif //USE_TENSORFLOW

	/// @brief Downloads model from GitLab if the specified path does not exist or is empty
	/// @param[in] path_to_model Directory path where the model should be located
	/// @param[in] base_url Url of the GitLab where models are stored
	/// @param[in] auto_download Whether to automatically download missing models
	static void download_model_if_not_existing( std::string const & path_to_model, std::string const & base_url, bool const & auto_download );

private: //Methods

	/// @brief Checks if the given model name is one of the available models
	static bool is_model_valid( std::string const & model_name) ;

	/// @brief Returns a string containing the names of all available models, separated by commas
	static std::string concatenate_model_names() ;

#ifdef USE_TENSORFLOW


#endif
private: //Data

	// defines the base directory where the models are store on the MeilerLab IWE GitLab
	static std::string const base_url_;
	// defines which models are available
	static std::vector< std::string > const available_models_;
	// defines whether models get automatically downloaded, default value of option is false
	bool auto_download_ = false;

protected:
#ifdef USE_TENSORFLOW
        /// @brief Set the tensorflow session.
        void set_tensorflow_session(basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session_in);
#endif // USE_TENSORFLOW
	/// @brief Given a selection, return a vector if selected indices.
	static utility::vector1< core::Size > get_selected_residues( core::select::residue_selector::ResidueSubset const & selected);


};

} //esm_perplexity
} //protocols


#endif //INCLUDED_protocols_esm_perplexity_EsmPerplexityTensorflowProtocol_hh

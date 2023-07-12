// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/esm_perplexity/PerResidueEsmProbabilitiesMetric.hh
/// @brief A class for predicting amino acid probabilities P(Sequence|Sequence) using the ESM language model.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adapted from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)

#ifndef INCLUDED_protocols_esm_perplexity_PerResidueEsmProbabilitiesMetric_HH
#define INCLUDED_protocols_esm_perplexity_PerResidueEsmProbabilitiesMetric_HH

#include <protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.fwd.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// citation manager
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// basic headers
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION
// C++ headers
#include <map>

namespace protocols {
namespace esm_perplexity {

///@brief A metric for estimating the amino acid probabilities using the ESM language model family
class PerResidueEsmProbabilitiesMetric : public core::simple_metrics::PerResidueProbabilitiesMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PerResidueEsmProbabilitiesMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PerResidueEsmProbabilitiesMetric( PerResidueEsmProbabilitiesMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PerResidueEsmProbabilitiesMetric() override;


public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///@brief Calculate the metric.
	std::map< core::Size, std::map< core::chemical::AA, core::Real >>
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Name of the class
	std::string
	name() const override;

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	///@brief Name of the metric
	std::string
	metric() const override;

public:

	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	core::simple_metrics::SimpleMetricOP
	clone() const override;

	/// @brief Get the ESM model that will be used for prediction.
	inline std::string model() const { return model_; }

	/// @brief Get the output filename of the pssm
	inline std::string write_pssm() const { return write_pssm_; }

	/// @brief Set the residue selector that we'll be using.
	/// @details Passing nullptr results in no residue selector being used.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Get the residue selector.
	/// @details If this returns nullptr, it means that no residue selector is being used.
	core::select::residue_selector::ResidueSelectorCOP residue_selector() const;

	/// @brief A second optional residue selector for attention masking
	void set_residue_selector2(core::select::residue_selector::ResidueSelectorCOP selector);

private: //Data

	/// @brief The ESM model to use for prediction.
	std::string model_;

	/// @brief The output filename of the pssm
	std::string write_pssm_;

	/// @brief Whether to multirun the network
	bool multirun_{};

	/// @brief An optional residue selector.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_ = nullptr;

	/// @brief The tensorflow protocol used.
	EsmPerplexityTensorflowProtocolCOP tensorflow_protocol_ = nullptr;

	/// @brief Residue selector for attention masking.
	core::select::residue_selector::ResidueSelectorCOP selector_two_ = nullptr;

	/// @brief Fill the return_map with the probabilities from the softmax_map for all amino acids
	/// @param[in] softmax_map The map containing all probabilities or logits returned from ESM (including special tokens)
	/// @param[in] return_map The map that will be returned by the calculate function, will be filled with probabilities/logits by this function
	static void fill_return_map(const std::map<core::Size, utility::vector1<core::Real>> &softmax_map,
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> &return_map);

	/// @brief Get the sequence of the selection
	std::string get_selection_sequence(std::string const &pose_sequence,
		std::map<core::Size, utility::vector1<core::Real>> const &position_map) const;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
	/// @brief This metric is unpublished.  It returns Moritz Ertelt as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

};

} //esm_perplexity
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_esm_perplexity_PerResidueEsmProbabilitiesMetric )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_esm_perplexity_PerResidueEsmProbabilitiesMetric_HH






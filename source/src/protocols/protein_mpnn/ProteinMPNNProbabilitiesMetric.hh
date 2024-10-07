// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_mpnn/ProteinMPNNProbabilitiesMetric.hh
/// @brief A metric for predicting amino acid probabilities using ProteinMPNN (Dauparas et al.).
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)


#ifndef INCLUDED_protocols_protein_mpnn_ProteinMPNNProbabilitiesMetric_HH
#define INCLUDED_protocols_protein_mpnn_ProteinMPNNProbabilitiesMetric_HH

// #include <protocols/protein_mpnn/EsmPerplexityTensorflowProtocol.fwd.hh>
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
namespace protein_mpnn {

///@brief A metric for estimating amino acid probabilities using the ProteinMPNN model
class ProteinMPNNProbabilitiesMetric : public core::simple_metrics::PerResidueProbabilitiesMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ProteinMPNNProbabilitiesMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ProteinMPNNProbabilitiesMetric( ProteinMPNNProbabilitiesMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ProteinMPNNProbabilitiesMetric() override;


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

	/// @brief Get the output filename of the pssm
	std::string write_pssm() const { return write_pssm_; }

	/// @brief Set the residue selector that we'll be using.
	/// @details Passing nullptr results in no residue selector being used.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Set residues whose coordinates should be passed to ProteinMPNN
	void
	set_coord_selector_rs( core::select::residue_selector::ResidueSelectorCOP coord_selector_rs );

	/// @brief Set residues whose sequence should be masked
	void
	set_sequence_mask_selector_rs( core::select::residue_selector::ResidueSelectorCOP sequence_mask_selector_rs );

	/// @brief Set tied positions. This is a collection of lists, each list containing residue selectors that should be tied.
	/// @details The first residues of each selector will be tied together, then the second, etc.
	///          Each residue selector must have the same number of residues.
	void
	set_tied_pos_rs( utility::vector1< utility::vector1< core::select::residue_selector::ResidueSelectorCOP > > tied_pos_rs );

	/// @brief Get the residue selector.
	/// @details If this returns nullptr, it means that no residue selector is being used.
	core::select::residue_selector::ResidueSelectorCOP residue_selector() const;

private: //Data

	/// @brief The output filename of the pssm
	std::string write_pssm_;

	/// @brief A vector of vectors of residue selectors to specify tied positions for ProteinMPNN
	utility::vector1< utility::vector1< core::select::residue_selector::ResidueSelectorCOP > > tied_pos_rs_;
	/// @brief A residue selector to specify backbone coordinates to be used as input for ProteinMPNN (default all)
	core::select::residue_selector::ResidueSelectorCOP coord_selector_rs_;
	/// @brief A residue selector to mask specific parts of the input sequence for ProteinMPNN
	core::select::residue_selector::ResidueSelectorCOP  sequence_mask_selector_rs_;

	/// @brief An optional residue selector to limit the returned probabilities to a selection
	core::select::residue_selector::ResidueSelectorCOP residue_selector_ = nullptr;

#ifdef USE_TORCH
	/// @brief Fill the return_map with the probabilities from the softmax_map for all amino acids
	/// @param[in] softmax_map The map containing all probabilities or logits
	/// @param[in] return_map The map that will be returned by the calculate function, will be filled with probabilities/logits by this function
	static void fill_return_map(const std::map<core::Size, utility::vector1<core::Real>> &softmax_map,
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> &return_map);
#endif // USE_TORCH
	/// @brief Get the sequence of the selection
	static std::string get_selection_sequence(std::string const &pose_sequence,
		std::map<core::Size, utility::vector1<core::Real>> const &position_map) ;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
	/// @brief This metric is currently unpublished, returns ME as author and the ProteinMPNN paper.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

};

} //protein_mpnn
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_protein_mpnn_ProteinMPNNProbabilitiesMetric )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_protein_mpnn_ProteinMPNNProbabilitiesMetric_HH






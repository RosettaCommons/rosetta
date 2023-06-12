// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ptm_prediction/PTMPredictionMetric.hh
/// @brief A class for predicting post-translational modifications using a Neural Net.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adapted from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)

#ifndef INCLUDED_protocols_ptm_prediction_PTMPredictionMetric_HH
#define INCLUDED_protocols_ptm_prediction_PTMPredictionMetric_HH

#include <protocols/ptm_prediction/PTMPredictionMetric.fwd.hh>
#include <protocols/ptm_prediction/PTMPredictionTensorflowProtocolBase.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

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
namespace ptm_prediction {

///@brief A metric for estimating the probability of a given residue(s) to be modified, as predicted by a pre-trained neural network.
class PTMPredictionMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PTMPredictionMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PTMPredictionMetric( PTMPredictionMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PTMPredictionMetric() override;


public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///@brief Calculate the metric.
	std::map< core::Size, core::Real >
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

	/// @brief Get the modification type that should be predicted
	inline std::string modification() const { return modification_; }

	/// @brief Set the residue selector that we'll be using.
	/// @details Passing nullptr results in no residue selector being used.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Get the residue selector.
	/// @details If this returns nullptr, it means that no residue selector is being used.
	core::select::residue_selector::ResidueSelectorCOP residue_selector() const;


	/// @brief enum class to define the different modifications that can be predicted
	enum class PTMPredictionMetricModes {
		INVALID_MODE = 0,
		ACETYLATION,
		ARG_METHYLATION,
		CITRULLINATION,
		CROTONYLATION,
		DEAMIDATION,
		GAMMA_CARBOXY_GLUTAMIC_ACID,
		GLUTARYLATION,
		GLUTATHIONYLATION,
		HYDROXYLATION,
		LYS_METHYLATION,
		MALONYLATION,
		N_LINKED_GLYCOSYLATION,
		O_LINKED_GLYCOSYLATION,
		PHOSPHORYLATION,
		S_NITROSYLATION,
		SUCCINYLATION,
		SUMOYLATION,
		UBIQUITINATION,
		N_MODES = UBIQUITINATION
	};

	/// @brief Get the modification type that should be predicted as enum class type
	inline PTMPredictionMetric::PTMPredictionMetricModes ptm_mode() const { return ptm_mode_; }

	/// @brief Convert an enum modification type into the respective string modification type
	static std::string enum_to_string( PTMPredictionMetric::PTMPredictionMetricModes & ptm_mode );

	/// @brief Convert an string modification type into an enum modification type
	static PTMPredictionMetric::PTMPredictionMetricModes string_to_enum( std::string & modification_type );

	// private:

private: //Data

	/// @brief The kind of post-translational modification to be predicted
	std::string modification_;

	/// @brief The kind of post-translational modification to be predicted as enum class
	PTMPredictionMetric::PTMPredictionMetricModes ptm_mode_;

	/// @brief An optional residue selector.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_ = nullptr;

	/// @brief The tensorflow protocol used to predict post-translational modifications.
	PTMPredictionTensorflowProtocolBaseCOP tensorflow_protocol_ = nullptr;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
	/// @brief This metric is unpublished.  It returns Moritz Ertelt as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;
};

} //ptm_prediction
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_ptm_prediction_PTMPredictionMetric )
#endif // SERIALIZATION


#endif //protocols_ptm_prediction_PTMPredictionMetric_HH


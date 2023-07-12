// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/esm_perplexity/PseudoPerplexityMetric.hh
/// @brief A class for calculating the pseudo-perplexity from a given PerResidueProbabilitiesMetric.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)
/// @note This has been adopted from how the ResidueSummaryMetric from Jared works.


#ifndef INCLUDED_protocols_esm_perplexity_EsmPerplexityMetric_HH
#define INCLUDED_protocols_esm_perplexity_EsmPerplexityMetric_HH

#include <protocols/esm_perplexity/PseudoPerplexityMetric.fwd.hh>
#include <protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>

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

///@brief A metric for calculating the (pseudo-)perplexity of a sequence using the ESM language model family
class PseudoPerplexityMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PseudoPerplexityMetric();

	PseudoPerplexityMetric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric );

	/// @brief Copy constructor (not needed unless you need deep copies)
	PseudoPerplexityMetric(PseudoPerplexityMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PseudoPerplexityMetric() override;


public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///@brief Calculate the metric.
	core::Real
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

	///@brief Set the PerResidueProbabilitiesMetric that will be used to calculate the pseudo-perplexity.
	void
	set_metric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric );

	///@brief Set a boolean to attempt to find cached data matching the name/custom_type of the passed in simple_metric.
	/// Optionally pass any set prefix/suffix.
	///
	/// This will allow the filter to re-use previously calculated data.
	///
	void
	set_use_cached_data( bool use_cache, std::string const & prefix="", std::string const & suffix="");

	///@brief If use_cache is set to false, do we fail if no data is found in the pose?
	/// Default True
	void
	set_fail_on_missing_cache(bool fail);

	/// @brief Function to return the (pseudo-)perplexity from a map of probabilities
	/// @details Takes a map of amino acid probabilities, gets the amino acid present in the pose, add the logarithm
	///          of its probability to a sum and then returns the exp of the negative ratio of this sum divided by the total length.
	/// @param[in] pose The pose used to refer which core::chemical::AA types are actually present
	/// @param[in] values The resulting values of a PerResidueProbabilitiesMetric
	/// @returns A core::Real value representing the pseud-perplexity
	static core::Real compute_perplexity(core::pose::Pose const &pose,
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> const &values) ;
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

	// metric to be used
	core::simple_metrics::PerResidueProbabilitiesMetricCOP metric_ = nullptr;

	//Accessing cached data from the pose
	bool use_cache_ = false;
	std::string cache_prefix_;
	std::string cache_suffix_;
	bool fail_on_missing_cache_ = true;

private: //Data


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
CEREAL_FORCE_DYNAMIC_INIT( protocols_esm_perplexity_EsmPerplexityMetric )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_esm_perplexity_EsmPerplexityMetric_HH






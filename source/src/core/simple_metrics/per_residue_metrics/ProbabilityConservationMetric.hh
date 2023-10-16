// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/ProbabilityConservationMetric.hh
/// @brief A class for calculating the conservation of a position given some predicted probabilities (using the relative Shannon Entropy).
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)


#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_ProbabilityConservationMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_ProbabilityConservationMetric_HH

#include <core/simple_metrics/per_residue_metrics/ProbabilityConservationMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

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

namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

///@brief A metric for calculating the conservation of a position given some predicted probabilities (using the relative Shannon Entropy).
class ProbabilityConservationMetric : public core::simple_metrics::PerResidueRealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ProbabilityConservationMetric();

	ProbabilityConservationMetric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric );

	/// @brief Copy constructor (not needed unless you need deep copies)
	ProbabilityConservationMetric(ProbabilityConservationMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ProbabilityConservationMetric() override;


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

	///@brief Set the PerResidueProbabilitiesMetric that will be used to calculate the conservation.
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

	///@brief calculate the delta
	std::map<core::Size, core::Real> compute_conservation( std::map<core::Size, std::map<core::chemical::AA, core::Real>> const & values  ) const;

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

} //per_residue_metrics
} //simple_metrics
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_ProbabilityConservationMetric )
#endif // SERIALIZATION


#endif //INCLUDED_core_simple_metrics_per_residue_metrics_ProbabilityConservationMetric_HH






// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/AverageProbabilitiesMetric.hh
/// @brief A metric for averaging multiple PerResidueProbabilitiesMetrics
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)


#ifndef INCLUDED_core_simple_metrics_metrics_AverageProbabilitiesMetric_HH
#define INCLUDED_core_simple_metrics_metrics_AverageProbabilitiesMetric_HH

#include <core/simple_metrics/metrics/AverageProbabilitiesMetric.fwd.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>

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
namespace metrics {

///@brief A metric for averaging multiple PerResidueProbabilitiesMetrics
class AverageProbabilitiesMetric : public core::simple_metrics::PerResidueProbabilitiesMetric {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	AverageProbabilitiesMetric();

	AverageProbabilitiesMetric(utility::vector1<core::simple_metrics::SimpleMetricCOP> metrics, utility::vector1< core::Real > weights);

	/// @brief Copy constructor (not needed unless you need deep copies)
	AverageProbabilitiesMetric(AverageProbabilitiesMetric const &src);

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~AverageProbabilitiesMetric() override;


public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///@brief Calculate the metric.
	std::map<core::Size, std::map<core::chemical::AA, core::Real>>
	calculate(core::pose::Pose const &pose) const override;

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

	///@brief Set the PerResidueProbabilitiesMetric that will be used to calculate.
	void
	set_metric(utility::vector1<core::simple_metrics::SimpleMetricCOP> metrics, utility::vector1< core::Real > weights );

	///@brief Set a boolean to attempt to find cached data matching the name/custom_type of the passed in simple_metric.
	/// Optionally pass any set prefix/suffix.
	///
	/// This will allow the filter to re-use previously calculated data.
	///
	void
	set_use_cached_data(bool use_cache, std::string const &prefix = "", std::string const &suffix = "");

	///@brief If use_cache is set to false, do we fail if no data is found in the pose?
	/// Default True
	void
	set_fail_on_missing_cache(bool fail);

public:

	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data) override;

	static
	void
	provide_xml_schema(utility::tag::XMLSchemaDefinition &xsd);

	core::simple_metrics::SimpleMetricOP
	clone() const override;

	///@brief computest the average between multiple PerResidueProbabilitiesMetrics outputs
	std::map<core::Size, std::map<core::chemical::AA, core::Real>>
	compute_average(utility::vector1<std::map<core::Size, std::map<core::chemical::AA, core::Real>>> const &all_values) const;

private: //Data
	// metrics to be averaged
	utility::vector1<core::simple_metrics::SimpleMetricCOP> metrics_;

	// weights factor for each metric
	utility::vector1< core::Real > weights_;

	//Accessing cached data from the pose
	bool use_cache_ = false;
	std::string cache_prefix_;
	std::string cache_suffix_;
	bool fail_on_missing_cache_ = true;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

	/// @brief This metric is unpublished.  It returns Moritz Ertelt as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList &citations) const override;

};

} //metrics
} //simple_metrics
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_metrics_AverageProbabilitiesMetric )
#endif // SERIALIZATION


#endif //INCLUDED_core_simple_metrics_metrics_AverageProbabilitiesMetric_HH






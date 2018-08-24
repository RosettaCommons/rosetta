// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/CompositeStringMetric.hh
///
/// @brief  Base class for core::Real CompositeStringMetrics
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_CompositeStringMetric_hh
#define INCLUDED_core_simple_metrics_CompositeStringMetric_hh

// Unit headers
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/CompositeStringMetric.fwd.hh>

// Core headers

#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ includes
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {

/// @brief A class that is derived to calculate a set of string values.
///  Apply(pose) method calculates this metric and adds it to the pose score for output.
///
/// @details
///  Calculate(pose) method calculates multiple string values and returns them as a map<string, string>
///   Name:Value
///
class CompositeStringMetric : public core::simple_metrics::SimpleMetric {

public: // constructors / destructors
	CompositeStringMetric();
	~CompositeStringMetric() override;

	CompositeStringMetric( CompositeStringMetric const & other );

	/// @brief Calculate the metric and add it to the pose within the SimpleMetricData cache.
	///           labeled as prefix+metric+suffix.
	///
	///
	/// @details  Data is output into the final score file, but can be accessed if needed through the cache.
	///
	void
	apply( pose::Pose & pose, std::string prefix="", std::string suffix="" ) const override;

	///@brief Calculate the metric.
	///
	///@details
	/// Return by value as this function can not STORE the result, it only calculates.
	/// Store the result in the pose by using the apply method, which calls this method and stores the result
	/// in the pose as ExtraScoreValues.
	virtual std::map< std::string, std::string >
	calculate( pose::Pose const & pose ) const = 0;

	///@brief Grab the data from the pose if it exists or calculate the metric
	///
	///@details If use_cache is true, we will attempt to pull the data from the pose.
	/// If fail_on_missing_cache is true, we will fail, otherwise, we will calculate the metric.
	///
	/// This function is meant to support caching metrics, so values do not need to be calculated twice,
	///  for example in SimpleMetricFilter/Features
	///  or code-wise where data takes a while to calculate and can be reused.
	///
	std::map< std::string, std::string >
	cached_calculate(
		pose::Pose const & pose,
		bool use_cache,
		std::string prefix="",
		std::string suffix="",
		bool fail_on_missing_cache=true) const;

public:

	///@brief Name of the class
	virtual std::string
	name() const override = 0;

	///@brief Name of the metric
	virtual std::string
	metric() const override = 0;

	///@brief Get the submetric names that this Metric will calculate
	// Used for features reporter.
	utility::vector1< std::string >
	get_metric_names() const override = 0;

public:
	/// @brief called by parse_my_tag -- should not be used directly
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override = 0;

	virtual SimpleMetricOP
	clone() const override = 0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class CompositeStringMetrics

} //namespace simple_metrics
} //namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_CompositeStringMetric )
#endif // SERIALIZATION

#endif

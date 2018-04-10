// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/CompositeRealMetric.hh
///
/// @brief  Base class for core::Real CompositeRealMetrics
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_CompositeRealMetric_hh
#define INCLUDED_core_simple_metrics_CompositeRealMetric_hh

// Unit headers
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/CompositeRealMetric.fwd.hh>

// Core headers

#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ includes
#include <map>

namespace core {
namespace simple_metrics {

/// @brief A class that is derived to calculate a set of core::Real values.
///  Apply(pose) method calculates this metric and adds it to the pose score for output.
///
///
/// @details
///  Calculate(pose) method calculates core::Real values and returns them as a map<std::string, Real>
///   Name:Value
///
class CompositeRealMetric : public core::simple_metrics::SimpleMetric {

public: // constructors / destructors

	CompositeRealMetric();

	~CompositeRealMetric() override;

	CompositeRealMetric( CompositeRealMetric const & other );

	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+calc_name+"_"+metric+suffix.
	///           calc_name is the individual component of the composite values calculated here
	///
	/// @details Score is added through setPoseExtraScore and is output
	///            into the score table/ score file at pose output.
	void
	apply( pose::Pose & pose, std::string prefix="", std::string suffix="" ) const override;

	///@brief Calculate the metric.
	virtual std::map< std::string, core::Real >
	calculate( pose::Pose const & pose ) const = 0;

public:

	///@brief Name of the class
	virtual std::string
	name() const override = 0;

	///@brief Name of the metric
	virtual std::string
	metric() const override = 0;

	///@brief Get the submetric names that this Metric will calculate

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

}; //class CompositeRealMetrics

} //namespace simple_metrics
} //namespace core


#endif

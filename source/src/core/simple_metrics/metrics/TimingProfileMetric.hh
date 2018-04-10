// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/TimingProfileMetric.hh
/// @brief Calculate the time difference between construction and apply/calculate.  Useful to time protocols in RosettaScripts or through mover containers.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_TimingProfileMetric_HH
#define INCLUDED_core_simple_metrics_metrics_TimingProfileMetric_HH

#include <core/simple_metrics/metrics/TimingProfileMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

//C++ headers
#include <chrono>

namespace core {
namespace simple_metrics {
namespace metrics {

///@brief Calculate the time difference between construction and apply/calculate.
/// Useful to time protocols in RosettaScripts or through mover containers.
///
class TimingProfileMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	TimingProfileMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	TimingProfileMetric( TimingProfileMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~TimingProfileMetric() override;

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///Defined in RealMetric:
	///
	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setExtraScorePose and is output
	///            into the score tables/file at pose output.
	//void
	//apply( pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Returns time from construction to the call of this function in minutes, with decimal.
	/// Options are available to calculate in hours instead.
	core::Real
	calculate( pose::Pose const & pose ) const override;

public:

	///@brief Set option to calculate the timings in hours.
	/// Default false (minutes)
	///
	void
	set_calc_in_hours( bool calc_in_hours );

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

private:
	std::chrono::high_resolution_clock::time_point construction_time_;
	bool calc_in_hours_ = false;
};

} //core
} //simple_metrics
} //metrics



#endif //core_simple_metrics_metrics_TimingProfileMetric_HH






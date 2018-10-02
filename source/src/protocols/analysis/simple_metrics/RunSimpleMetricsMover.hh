// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/RunSimpleMetricsMover.hh
/// @brief Run the set of SimpleMetrics.  Add the data to the pose, which will be output int the score file.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_analysis_simple_metrics_RunSimpleMetricsMover_hh
#define INCLUDED_protocols_analysis_simple_metrics_RunSimpleMetricsMover_hh

// Unit headers
#include <protocols/analysis/simple_metrics/RunSimpleMetricsMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace analysis {
namespace simple_metrics {

///@brief Run the set of SimpleMetrics.  Add the data to the pose, which will be output int the score file.
///  with the given prefix and suffix. (prefix+metric_name+suffix)
///
class RunSimpleMetricsMover : public protocols::moves::Mover {

public:
	RunSimpleMetricsMover();
	RunSimpleMetricsMover( utility::vector1< core::simple_metrics::SimpleMetricCOP> const & metrics );

	RunSimpleMetricsMover( RunSimpleMetricsMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	~RunSimpleMetricsMover() override;

	///@brief Add a simple metric to run.
	void
	add_simple_metric( core::simple_metrics::SimpleMetricCOP metric );

	void
	set_simple_metrics( utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics );

	///@brief Run the set of SimpleMetrics.  Add the data to the pose, which will be output int the score file.
	///  with the given prefix and suffix. (prefix+metric_name+suffix)
	///
	void
	apply( core::pose::Pose & pose, std::string const & prefix, std::string const & suffix);

	///@brief Run the set of SimpleMetrics.  Add the data to the Pose, which will be output in the score file.
	void
	apply( core::pose::Pose & pose ) override;

public:

	void
	set_prefix( std::string const & prefix );

	void
	set_suffix( std::string const & suffix );

	///@brief Should we override existing data?
	/// Default False.
	///
	///@details
	///    If existing data is found and this option is false, we will throw an exception.
	///
	void
	set_override( bool override_existing_data);

public:

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

public:

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics_; //They do not accumulate state.
	std::string prefix_="";
	std::string suffix_="";
	bool override_existing_data_ = false;

};

} //simple_metrics
} //analysis
} //protocols

#endif //protocols_analysis_simple_metrics_RunSimpleMetricsMover_hh

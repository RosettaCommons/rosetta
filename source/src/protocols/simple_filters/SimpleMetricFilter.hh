// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SimpleMetricFilter.hh
/// @brief A filter takes any RealMetric and applies a set cutoff to filter the model.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_SimpleMetricFilter_hh
#define INCLUDED_protocols_simple_filters_SimpleMetricFilter_hh

// Unit headers
#include <protocols/simple_filters/SimpleMetricFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/simple_metrics/RealMetric.fwd.hh>
#include <core/simple_metrics/SimpleMetric.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

#include <string>
#include <map>

//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace simple_filters {

///@brief Enum that tells us how to interpret the cutoff set.
enum comparison_type{
	eq = 1,
	ne = 2,
	lt = 3,
	gt = 4,
	lt_or_eq = 5,
	gt_or_eq = 6,
	bogus = 7,
	comparison_type_total = bogus

};

/// @brief A filter takes any SimpleMetric and applies a set cutoff to filter the model.
///  Set the cutoff type to set the behavior of the metric.
///
/// @details
///  Compares metric_value to cutoff_value or match_value
///
///  match_value is used for strings, cutoff is used for numbers.
///
///  So, if comparison_type is set to eq, if metric_value eq the cutoff_value we return true.
///
///  RMSD example: set to lt, and the filter will pass with a metric value anything less than the cuttoff.
///
/// CompositeMetrics take an extra parameter, composite_action.
///  Composite action can be "any", "all", or any specific composite metric value you want to match on.
///
///  Any:  If any match the set criteria, we return true
///  All:  If all match the set criteria, we return true
///  etc:  If the particular metric of the composite metric (ex fa_rep for the CompositeEnergyMetric), matches the criteria,
///          we return true.
///
/// See Methods:
///  set_comparison_type(), set_composite_action(),
///  set_cutoff_value(), set_match_string()
///
class SimpleMetricFilter : public protocols::filters::Filter {

public:
	SimpleMetricFilter();

	SimpleMetricFilter(
		core::simple_metrics::SimpleMetricCOP metric,
		comparison_type co_type);


	SimpleMetricFilter(
		core::simple_metrics::SimpleMetricCOP metric,
		comparison_type co_type,
		std::string composite_action);


	// destructor (important for properly forward-declaring smart-pointer members)
	~SimpleMetricFilter() override;

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

public:

	///@brief Set the SimpleMetric that we will be using to filter.
	void
	set_simple_metric( core::simple_metrics::SimpleMetricCOP metric );

public:

	///@brief Set the cutoff value for any RealMetric or CompositeRealMetric.
	void
	set_cutoff( core::Real cutoff );


	///@brief Set the match value for any StringMetric or CompositeStringMetric.
	void
	set_match_string( std::string const & match_string );

public:

	///@brief Sets the cutoff type - aka eq, ne, etc. Options are:
	///  eq, ne, lt, gt, lt_or_eq, gt_or_eq.
	///
	///@details If this is a StringMetric, only eq and ne are relevant here.
	/// IF value [comparison_type] cutoff_ return True.
	void
	set_comparison_type( comparison_type co_type );

	///@brief
	///
	///  Set the action we take for a set CompositeMetric
	///
	///  Composite action can be "any", "all", or any specific composite metric value you want to match on.
	///  See Also: set_comparison_type(), set_cutoff
	///
	///@details
	///
	///  Any:  If any match the set criteria, we return true
	///  All:  If all match the set criteria, we return true
	///  etc:  If the particular metric of the composite metric (ex fa_rep for the CompositeEnergyMetric), matches the criteria,
	///          we return true.
	///
	void
	set_composite_action( std::string const & composite_action );

	///@brief Set the filter to use the SUM of values from a PerResidueRealMetric for filtering.
	///  This is instead of using the metric as a Composite Metric for each resnum.
	/// Default False.
	///
	void
	set_use_sum_for_per_residue_real( bool use_sum_for_per_residue_real );

public:

	///@brief Set the sigfigs we will use in our comparisons.
	/// Default is .0001;
	void
	set_epsilon( core::Real epsilon );

public:
	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

private:

	///Takes std::map of std::string - value and returns a boolean if it passes the filter.
	template <class T>
	bool compare_composites( std::map< std::string, T > const & values) const;

	///Takes std::map of std::string - value and returns a boolean if it passes the filter.
	template <class T>
	bool compare_composites( std::map< core::Size, T > const & values) const;

	bool
	compare_metric( core::Real value ) const ;

	bool
	compare_metric( std::string const & value ) const ;

private:

	core::simple_metrics::SimpleMetricCOP metric_ = nullptr;
	comparison_type type_ = bogus;
	core::Real cutoff_ = 0;
	std::string match_ = "";

	core::Real epsilon_ = .0001;
	std::string composite_action_ = "";
	bool sum_per_residue_real_metric_ = false;

};

std::map< std::string, comparison_type >
get_string_comparison_type_map();

utility::vector1< std::string >
get_string_comparison_type_strings();

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_SimpleMetricFilter_hh

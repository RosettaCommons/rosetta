// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/ResidueSummaryMetric.hh
/// @brief A SimpleMetric that takes a PerResidueMetric and calculates different summaries of the overall data.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_ResidueSummaryMetric_HH
#define INCLUDED_core_simple_metrics_metrics_ResidueSummaryMetric_HH

#include <core/simple_metrics/metrics/ResidueSummaryMetric.fwd.hh>
#include <core/simple_metrics/PerResidueRealMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace simple_metrics {
namespace metrics {

///@brief Enum that tells us HOW to summarize the PerResidueMetric data
enum summary_type{
	sum = 1,
	mean,
	n_res_eq,
	n_res_ne,
	n_res_gt,
	n_res_lt,
	n_res_gt_or_eq,
	n_res_lt_or_eq,
	bogus,

	summary_type_total = bogus

};


///@brief A SimpleMetric that takes a PerResidueRealMetric and calculates different summaries of the overall data.
///
///@details
/// Can compute totals, means, or the number of residues matching a certain criteria,
//// for example less than or equal to a certain value.
///  Useful for summarizing metrics or using more complex functions for filters.
///
class ResidueSummaryMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ResidueSummaryMetric();

	ResidueSummaryMetric( PerResidueRealMetricCOP metric );

	/// @brief Copy constructor (not needed unless you need deep copies)
	ResidueSummaryMetric( ResidueSummaryMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ResidueSummaryMetric() override;

public:

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
	//apply( core::pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Calculate the metric.
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Set the SimpleMetric that we will be using to filter.
	void
	set_metric( core::simple_metrics::PerResidueRealMetricCOP metric );

	///@brief Set the summary behavior of this class.
	/// Default mean
	void
	set_action( summary_type sum_type);

	///@brief Set any cutoff for summing n_residues
	/// Default 0
	void
	set_action_value( core::Real cutoff );

public:

	///@brief Set a boolean to attempt to find cached data matching the name/custom_type of the passed in simple_metric.
	/// Optionally pass any set prefix/suffix.
	///
	/// This will allow the filter to re-use previously calculated data.
	///
	void
	set_use_cached_data( bool use_cache, std::string prefix="", std::string suffix="");

	///@brief If use_cache is set to false, do we fail if no data is found in the pose?
	/// Default True
	void
	set_fail_on_missing_cache(bool fail);




public:

	///@brief Set the sigfigs we will use in our comparisons.
	/// Default is .0001;
	///
	///@details
	///  The absolute value of the difference between two numbers, below which they will be considered to be equal.
	///   Used when comparing numbers if action = n_res_eq or n_res_ne.
	void
	set_epsilon( core::Real epsilon );

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

private:

	PerResidueRealMetricCOP metric_ = nullptr;

	core::Real cutoff_ = 0;
	core::Real epsilon_ = .0001;
	summary_type action_ = mean;

	//Accessing cached data from the pose
	bool use_cache_ = false;
	std::string cache_prefix_ = "";
	std::string cache_suffix_ = "";
	bool fail_on_missing_cache_ = true;

};

static const std::map< std::string, summary_type > summary_string_to_type(
	{
	{ "sum", sum },
	{ "mean", mean },
	{ "n_res_eq", n_res_eq },
	{ "n_res_ne", n_res_ne },
	{ "n_res_lt", n_res_lt },
	{ "n_res_gt", n_res_gt },
	{ "n_res_lt_or_eq", n_res_lt_or_eq },
	{ "n_res_gt_or_eq", n_res_gt_or_eq }
	}
);

utility::vector1< std::string >
get_summary_type_strings();


} //core
} //simple_metrics
} //metrics

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_metrics_ResidueSummaryMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_metrics_ResidueSummaryMetric_HH






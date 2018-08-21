// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/composite_metrics/ProtocolSettingsMetric.hh
/// @brief This Metric reports options that have been set in the command line and splits script_vars.  Each option name is the type and the setting is the value in the map.  This is primarily aimed at benchmarking and record-keeping for large-scale rosetta runs or experiments.  It works with both the global and local OptionsCollection to enable its use in JD3.  It is analogous to the ProtocolFeatures reporter, with more options for xml-based variables.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_composite_metrics_ProtocolSettingsMetric_HH
#define INCLUDED_core_simple_metrics_composite_metrics_ProtocolSettingsMetric_HH



#include <core/simple_metrics/composite_metrics/ProtocolSettingsMetric.fwd.hh>
#include <core/simple_metrics/CompositeStringMetric.hh>

// Core headers
#include <core/types.hh>


// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/options/OptionCollection.fwd.hh>

// C++ headers
#include <map>

namespace core {
namespace simple_metrics {
namespace composite_metrics {

///@brief This Metric reports options that have been set in the command line and splits script_vars.  Each option name is the type and the setting is the value in the map.  This is primarily aimed at benchmarking and record-keeping for large-scale rosetta runs or experiments.  It works with both the global and local OptionsCollection to enable its use in JD3.  It is analogous to the ProtocolFeatures reporter, with more options for xml-based variables.
class ProtocolSettingsMetric : public core::simple_metrics::CompositeStringMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Parse the Global Options Collections
	ProtocolSettingsMetric(
		bool base_name_option_only = true,
		bool get_script_vars = true,
		bool get_user_options = true,
		bool skip_corrections = true);

	///@brief Parse the Local Options Collection
	ProtocolSettingsMetric(
		utility::options::OptionCollection const & options,
		bool base_name_option_only = true,
		bool get_script_vars = true,
		bool get_user_options = true,
		bool skip_corrections = true );


	/// @brief Copy constructor (not needed unless you need deep copies)
	ProtocolSettingsMetric( ProtocolSettingsMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ProtocolSettingsMetric() override;

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

	///@brief Calculate the metric.
	std::map< std::string, std::string >
	calculate( core::pose::Pose const & pose ) const override;

	///@brief Only return these options.  If not getting full namespaces, these are ONLY the base names.
	void
	set_only_report_these_options(utility::vector1< std::string > const & select_opts);

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

	///@brief Get the submetric names that this Metric will calculate
	utility::vector1< std::string >
	get_metric_names() const override;

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

public:

	///@brief Parse an options collection into a set of string-string pairs.
	void
	parse_options(
		utility::options::OptionCollection const & options,
		bool base_name_option_only = true,
		bool get_script_vars = true,
		bool get_user_options = true,
		bool skip_corrections = true);

	///@brief Split the script_vars string into individual options and place them into options_values map.
	void
	split_script_vars(std::string const & script_vars_option_string, std::map< std::string, std::string> & options_values) const;

private:
	std::map< std::string, std::string > options_values_;

	utility::vector1< std::string > limit_reporting_to_these_options_;
};


} //core
} //simple_metrics
} //composite_metrics



#endif //core_simple_metrics_composite_metrics_ProtocolSettingsMetric_HH






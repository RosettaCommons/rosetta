// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/CalculatorMetric.hh
/// @brief A metric which can combine other metrics in a (semi) arbitrary calculation
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_CalculatorMetric_HH
#define INCLUDED_core_simple_metrics_metrics_CalculatorMetric_HH

#include <core/simple_metrics/metrics/CalculatorMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <numeric/Calculator.fwd.hh>
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {
namespace metrics {

/// @brief A metric which can combine other metrics
class CalculatorMetric : public core::simple_metrics::RealMetric {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	CalculatorMetric();

	/// @brief Initialize with equation.
	CalculatorMetric(std::string const & equation);

	/// @brief For the particular variable name in the equation, use the passed metric.
	/// The passed metric must be a RealMetric, else you'll get an error message.
	void add_simple_metric( std::string const & name, core::simple_metrics::SimpleMetricCOP metric );

	/// @brief For the variable name in the equation, look up the report_key scoreterm stored in the Pose
	void add_reported_value( std::string const & name, std::string const & report_key );

	/// @brief For the varible name in the equation, use the passed constant value
	void add_constant( std::string const & name, core::Real value );

	/// @brief Do a quick check to make sure the equation parsing worked out.
	/// Returns true on success, false on failure
	bool check_equation();

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

	/// @brief Calculate the metric.
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

public:

	/// @brief Name of the class
	std::string
	name() const override;

	/// @brief Name of the class for creator.
	static
	std::string
	name_static();

	/// @brief Name of the metric
	std::string
	metric() const override;

	/// @brief This simple metric is unpublished (returns true).
	bool simple_metric_is_unpublished() const override;

	/// @brief This simple metric is unpublished.
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
	provide_authorship_info_for_unpublished() const override;

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
	numeric::CalculatorOP calc_;
	std::map<std::string, core::Real> values_;
	std::map<std::string, core::simple_metrics::RealMetricCOP> metrics_;
	std::map<std::string, std::string> reported_values_;

};

} //metrics
} //simple_metrics
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_metrics_CalculatorMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_metrics_CalculatorMetric_HH






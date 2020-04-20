// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/FilterValueMetric.hh
/// @brief Convert the result of a Filter's report_sm() to a SimpleMetric
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_filters_FilterValueMetric_HH
#define INCLUDED_protocols_filters_FilterValueMetric_HH

#include <protocols/filters/FilterValueMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace filters {

/// @brief Convert the result of a Filter's report_sm() to a SimpleMetric
/// This is intended primarily as a compatibility shim class for making old code compatible.
/// It's not really intended to be used for new work - write a SimpleMetric directly.
class FilterValueMetric : public core::simple_metrics::RealMetric {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	FilterValueMetric();

	/// @brief Construct based on a given filter.
	FilterValueMetric( FilterOP filter );

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

	protocols::filters::FilterOP filter_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //filters
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_filters_FilterValueMetric )
#endif // SERIALIZATION

#endif //protocols_filters_FilterValueMetric_HH






// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/CustomStringValueMetric.hh
/// @brief A simple metric that allows an arbitrary, user- or developer-set string to be cached in a pose.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_core_simple_metrics_metrics_CustomStringValueMetric_HH
#define INCLUDED_core_simple_metrics_metrics_CustomStringValueMetric_HH

#include <core/simple_metrics/metrics/CustomStringValueMetric.fwd.hh>
#include <core/simple_metrics/StringMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {
namespace metrics {

///@brief A simple metric that allows an arbitrary, user- or developer-set string to be cached in a pose.
class CustomStringValueMetric : public core::simple_metrics::StringMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	CustomStringValueMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	CustomStringValueMetric( CustomStringValueMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~CustomStringValueMetric() override;


public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///Defined in StringMetric:
	///
	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setExtraScorePose and is output
	///            into the score tables/file at pose output.
	//void
	//apply( core::pose::Pose & pose, prefix="", suffix="" ) override;

	/// @brief Calculate the metric.
	/// @details Returns the cached value.
	std::string
	calculate( core::pose::Pose const & pose ) const override;

	/// @brief Set the value that we're going to cache in the pose.
	inline void set_value( std::string const &value_in ) { value_ = value_in; }

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

private: //Data

	/// @brief The cached value.
	std::string value_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //core
} //simple_metrics
} //metrics

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_metrics_CustomStringValueMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_metrics_CustomStringValueMetric_HH






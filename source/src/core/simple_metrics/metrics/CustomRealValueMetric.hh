// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/CustomRealValueMetric.hh
/// @brief A simple metric that allows an arbitrary, user- or developer-set floating-point value to be cached in a pose.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_core_simple_metrics_metrics_CustomRealValueMetric_HH
#define INCLUDED_core_simple_metrics_metrics_CustomRealValueMetric_HH

#include <core/simple_metrics/metrics/CustomRealValueMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

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

///@brief A simple metric that allows an arbitrary, user- or developer-set floating-point value to be cached in a pose.
class CustomRealValueMetric : public core::simple_metrics::RealMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	CustomRealValueMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	CustomRealValueMetric( CustomRealValueMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~CustomRealValueMetric() override;

public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	/// @brief Calculate the metric.
	/// @details Returns the cached value.
	core::Real
	calculate(
		core::pose::Pose const & pose
	) const override;

	/// @brief Set the value that we're going to cache in the pose.
	inline void set_value( core::Real const value_in ) { value_ = value_in; }

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
		basic::datacache::DataMap & data
	) override;

	static
	void
	provide_xml_schema(
		utility::tag::XMLSchemaDefinition & xsd
	);

	core::simple_metrics::SimpleMetricOP
	clone() const override;

public: //Functions needed for the citation manager

	/// @brief Provide the citation.
	void
	provide_citation_info(
		basic::citation_manager::CitationCollectionList & citations
	) const override;

private: //Data

	/// @brief The cached value.
	core::Real value_ = 0.0;

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
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_metrics_CustomRealValueMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_metrics_CustomRealValueMetric_HH






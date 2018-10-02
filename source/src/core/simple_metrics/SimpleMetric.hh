// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/AbstractMetric.hh
/// @brief The base class for Metrics in the Metric/Filter/Reporter system
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @note   This is interface: it has no fields, and only
///         pure virtual methods.  No further constructors should
///         be defined.


#ifndef INCLUDED_core_simple_metrics_SimpleMetric_hh
#define INCLUDED_core_simple_metrics_SimpleMetric_hh


// Project forward headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {


/// @brief The base class for Metrics in the Metric/Filter/Reporter system
/// @details The non-templated base class allows us to build one from a factory and interact with it
/// through RosettaScripts
class SimpleMetric : public utility::pointer::ReferenceCount {

public:


	SimpleMetric( std::string const & simple_metric_type );

	virtual
	~SimpleMetric();

	SimpleMetric( SimpleMetric const & other );

	//Every SimpleMetric should implement a calculate method.
	// Due to the use of a factory and owning pointers, this is not added to the abstract base class.  But please, be aware.
	// Currently, this code here does nothing useful.
	//template < typename T >
	//T calculate( const pose::Pose & pose ) const;

	///@brief Calculate the metric and add it to the Score, which is output into a scorefile - labeled as prefix+metric+suffix.
	/// Must be implemented by derived classes
	virtual void
	apply(
		pose::Pose & pose,
		std::string prefix="",
		std::string suffix="",
		bool override_existing_data=false) const = 0;

	///@brief Get the name of SimpleMetric class
	virtual std::string
	name() const = 0;

	///@brief Get the name of the Metric
	virtual std::string
	metric() const = 0;

	virtual SimpleMetricOP
	clone() const = 0;

	///@brief Get the metric name(s) that this Metric will calculate
	virtual utility::vector1< std::string >
	get_metric_names() const = 0;

	void
	set_custom_type( std::string const & custom_type );

	///@brief Additional setting to prefix/suffix
	//  so that many different configured SMs can be called in one RunSimpleMetric run
	/// Output data name will be prefix+custom_type+type+suffix
	std::string
	get_custom_type() const;

public:

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) = 0;

	///Parse the base class tag.  Keep required interface for parse_my_tag.
	virtual void
	parse_base_tag(
		utility::tag::TagCOP tag );

	static
	utility::tag::XMLSchemaComplexTypeGeneratorOP
	complex_type_generator_for_simple_metric( utility::tag::XMLSchemaDefinition & );

	std::string
	simple_metric_type() const {
		return simple_metric_type_;
	};

	///@brief Get the final name of this metric including its simple_metric_type_ name and any set custom type.
	///
	std::string
	get_final_sm_type() const;

private:

	///@brief Type of SimpleMetric.  AKA RealMetric, StringMetric, etc.
	std::string simple_metric_type_;

	//Additional setting to prefix/suffix for RosettaScripts -
	//  so that many different configured SMs can be called in one RunSimpleMetric run
	std::string custom_type_ = "";

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // SimpleMetric


} //core
} //simple_metrics

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_SimpleMetric )
#endif // SERIALIZATION


#endif //INCLUDED_core_metrics_AbstractMetric_hh




// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/CustomStringValueMetric.cc
/// @brief A simple metric that allows an arbitrary, user- or developer-set string to be cached in a pose.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <core/simple_metrics/metrics/CustomStringValueMetric.hh>
#include <core/simple_metrics/metrics/CustomStringValueMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.metrics.CustomStringValueMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CustomStringValueMetric::CustomStringValueMetric():
	core::simple_metrics::StringMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
CustomStringValueMetric::CustomStringValueMetric( CustomStringValueMetric const &  ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CustomStringValueMetric::~CustomStringValueMetric(){}

core::simple_metrics::SimpleMetricOP
CustomStringValueMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new CustomStringValueMetric( *this ) );

}

std::string
CustomStringValueMetric::name() const {
	return name_static();
}

std::string
CustomStringValueMetric::name_static() {
	return "CustomStringValueMetric";

}
std::string
CustomStringValueMetric::metric() const {
	return get_custom_type().empty() ? "custom_string_valued_metric" : get_custom_type();
}

void
CustomStringValueMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  )
{
	SimpleMetric::parse_base_tag( tag );
	set_value( tag->getOption< std::string >( "value", "" ) );
}

void
CustomStringValueMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("value", xs_string, "The custom string to cache in the pose.", "");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for caching an arbitrary, user-defined string in a pose.", attlist);
}

std::string
CustomStringValueMetric::calculate(const core::pose::Pose & ) const {
	return value_;
}

void
CustomStringValueMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	CustomStringValueMetric::provide_xml_schema( xsd );
}

std::string
CustomStringValueMetricCreator::keyname() const {
	return CustomStringValueMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
CustomStringValueMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new CustomStringValueMetric );
}

} //core
} //simple_metrics
} //metrics


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::metrics::CustomStringValueMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::StringMetric>( this ) );
	arc( CEREAL_NVP( value_ ) );

}

template< class Archive >
void
core::simple_metrics::metrics::CustomStringValueMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::StringMetric >( this ) );
	arc( value_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::metrics::CustomStringValueMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::metrics::CustomStringValueMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_metrics_CustomStringValueMetric )
#endif // SERIALIZATION





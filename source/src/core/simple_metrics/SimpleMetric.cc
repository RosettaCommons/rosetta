// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/simple_metrics/SimpleMetric.cc
/// @brief The base class for Metrics in the Metric/Filter/Reporter system
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Unit header or inline function header
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/util.hh>

// NOTE: This file should have NO dependencies other than its header.

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace simple_metrics {



SimpleMetric::SimpleMetric( std::string const & simple_metric_type_name ):
	utility::pointer::ReferenceCount(),
	simple_metric_type_( simple_metric_type_name )

{

}

SimpleMetric::~SimpleMetric(){}

SimpleMetric::SimpleMetric( SimpleMetric const & src ):
	utility::pointer::ReferenceCount(),
	simple_metric_type_(src.simple_metric_type_)
{

}

void
SimpleMetric::set_custom_type(std::string const & custom_type){
	custom_type_ = custom_type;
}

std::string
SimpleMetric::get_custom_type() const {
	return custom_type_;
}

void
SimpleMetric::parse_base_tag(utility::tag::TagCOP tag ){
	set_custom_type(tag->getOption< std::string >("custom_type", custom_type_));
}

std::string
SimpleMetric::get_final_sm_type() const{
	std::string custom_type = get_custom_type();

	if ( custom_type != "" ) custom_type=custom_type+"_";
	return custom_type + metric();
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
SimpleMetric::complex_type_generator_for_simple_metric( utility::tag::XMLSchemaDefinition &  ) {

	using namespace utility::tag;

	AttributeList attlist;

	std::string custom_type_descrition =
		"Allows multiple configured SimpleMetrics of a single type to be called in a single RunSimpleMetrics and SimpleMetricFeatures."
		" \n The custom_type name will be added to the data tag in the scorefile or features database.";

	attlist
		+ XMLSchemaAttribute( "custom_type", xs_string, custom_type_descrition );


	XMLSchemaComplexTypeGeneratorOP ct_gen( new utility::tag::XMLSchemaComplexTypeGenerator );
	ct_gen->
		add_attributes( attlist )
		.complex_type_naming_func( & complex_type_name_for_simple_metric );


	return ct_gen;

}




} //core
} //metrics

#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::SimpleMetric::save( Archive & arc ) const {
	arc( CEREAL_NVP( simple_metric_type_));
	arc( CEREAL_NVP( custom_type_ ) );

}

template< class Archive >
void
core::simple_metrics::SimpleMetric::load( Archive & arc ) {
	arc( simple_metric_type_ );
	arc( CEREAL_NVP( custom_type_ ));
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::SimpleMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::SimpleMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_SimpleMetric )
#endif // SERIALIZATION

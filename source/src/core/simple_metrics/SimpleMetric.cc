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

utility::tag::XMLSchemaComplexTypeGeneratorOP
SimpleMetric::complex_type_generator_for_simple_metric( utility::tag::XMLSchemaDefinition &  ) {

	using namespace utility::tag;

	AttributeList attlist;

	std::string custom_type_descrition = "Additional setting to prefix/suffix "
		"so that many different configured SMs can be called in one RunSimpleMetric run\n"
		"  Output data name will be prefix+custom_type+type+suffix";

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



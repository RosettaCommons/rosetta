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

#include <basic/options/keys/simple_metrics.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// NOTE: This file should have NO dependencies other than its header.

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace simple_metrics {
using namespace basic::options;


SimpleMetric::SimpleMetric( std::string const & simple_metric_type_name ):
	utility::VirtualBase(),
	simple_metric_type_( simple_metric_type_name )

{

}

SimpleMetric::~SimpleMetric(){}

SimpleMetric::SimpleMetric( SimpleMetric const & src ):
	utility::VirtualBase(),
	simple_metric_type_(src.simple_metric_type_),
	custom_type_(src.custom_type_)
{}

SimpleMetric &
SimpleMetric::operator=( SimpleMetric const & ot ) {
	utility::VirtualBase::operator=( ot );
	simple_metric_type_ = ot.simple_metric_type_;
	custom_type_ = ot.custom_type_;
	return *this;
}

void
SimpleMetric::apply(
	pose::Pose & pose,
	std::string const & prefix,
	std::string const & suffix,
	bool override_existing_data
) const {
	apply( prefix + get_final_sm_type() + suffix, pose, override_existing_data );
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

	if ( tag->hasOption("name") ) {
		tag_name_ = tag->getOption< std::string >("name");
	}

	if ( option [OptionKeys::simple_metrics::use_name_as_custom_type]() && ! tag_name_.empty() ) {
		set_custom_type(tag_name_);
	}
	set_custom_type(tag->getOption< std::string >("custom_type", custom_type_));
}

std::string
SimpleMetric::get_final_sm_type() const{
	std::string custom_type = get_custom_type();

	if ( custom_type != "" ) custom_type=custom_type+"_";

	if ( option [OptionKeys::simple_metrics::use_name_as_metric]() && ! tag_name_.empty() ) {
		return tag_name_;
	} else {
		return custom_type + metric();
	}
}

/// @brief Provide citations to the passed CitationCollectionList
/// Subclasses should add the info for themselves and any other classes they use.
/// @details The default implementation of this function does nothing.  It may be
/// overriden by simple metrics wishing to provide citation information.
void
SimpleMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & ) const {
	// Do nothing
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

void
run_metrics( core::pose::Pose & pose, utility::vector1< SimpleMetricCOP > const & metrics, std::string prefix, std::string suffix){
	for ( SimpleMetricCOP const & metric : metrics ) {
		metric->apply( pose, prefix, suffix );
	}
}


} //core
} //metrics

#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::SimpleMetric::save( Archive & arc ) const {
	arc( CEREAL_NVP( simple_metric_type_));
	arc( CEREAL_NVP( custom_type_ ) );
	arc( CEREAL_NVP( tag_name_ ) );
}

template< class Archive >
void
core::simple_metrics::SimpleMetric::load( Archive & arc ) {
	arc( simple_metric_type_ );
	arc( custom_type_ );
	arc( tag_name_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::SimpleMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::SimpleMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_SimpleMetric )
#endif // SERIALIZATION

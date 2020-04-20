// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/FilterValueMetric.cc
/// @brief Convert the result of a Filter's report_sm() to a SimpleMetric
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/filters/FilterValueMetric.hh>
#include <protocols/filters/FilterValueMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.filters.FilterValueMetric" );


namespace protocols {
namespace filters {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
FilterValueMetric::FilterValueMetric():
	core::simple_metrics::RealMetric()
{}

FilterValueMetric::FilterValueMetric( FilterOP filter ):
	core::simple_metrics::RealMetric(),
	filter_( filter )
{}

core::simple_metrics::SimpleMetricOP
FilterValueMetric::clone() const {
	return utility::pointer::make_shared< FilterValueMetric >( *this );
}

std::string
FilterValueMetric::name() const {
	return name_static();
}

std::string
FilterValueMetric::name_static() {
	return "FilterValueMetric";

}
std::string
FilterValueMetric::metric() const {
	return "foo";
}

void
FilterValueMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data)
{
	SimpleMetric::parse_base_tag( tag );

	filter_ = protocols::rosetta_scripts::parse_filter( "filter", data );
}

void
FilterValueMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("filter", xs_string, "The filter to convert into a SimpleMetric.");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric which converts the result of a Filter::report_sm() calculation into a RealMetric. "
		"This is intended primarily as a compatibility shim, and not for greenfield code."
		, attlist);
}

core::Real
FilterValueMetric::calculate(const core::pose::Pose & pose ) const {

	if ( filter_ == nullptr ) {
		utility_exit_with_message("Filter to use in FilterValueMetric not set.");
	} else {
		return filter_->report_sm( pose );
	}
}



void
FilterValueMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	FilterValueMetric::provide_xml_schema( xsd );
}

std::string
FilterValueMetricCreator::keyname() const {
	return FilterValueMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
FilterValueMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< FilterValueMetric >();
}

} //filters
} //protocols


#ifdef    SERIALIZATION



template< class Archive >
void
protocols::filters::FilterValueMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	arc( CEREAL_NVP( filter_ ) );

}

template< class Archive >
void
protocols::filters::FilterValueMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );
	arc( filter_ );


}

SAVE_AND_LOAD_SERIALIZABLE( protocols::filters::FilterValueMetric );
CEREAL_REGISTER_TYPE( protocols::filters::FilterValueMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_filters_FilterValueMetric )
#endif // SERIALIZATION





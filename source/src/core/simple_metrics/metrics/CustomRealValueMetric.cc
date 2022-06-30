// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/CustomRealValueMetric.cc
/// @brief A simple metric that allows an arbitrary, user- or developer-set floating-point value to be cached in a pose.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <core/simple_metrics/metrics/CustomRealValueMetric.hh>
#include <core/simple_metrics/metrics/CustomRealValueMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.metrics.CustomRealValueMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CustomRealValueMetric::CustomRealValueMetric():
	core::simple_metrics::RealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
CustomRealValueMetric::CustomRealValueMetric( CustomRealValueMetric const &  ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CustomRealValueMetric::~CustomRealValueMetric(){}

core::simple_metrics::SimpleMetricOP
CustomRealValueMetric::clone() const {
	return utility::pointer::make_shared< CustomRealValueMetric >( *this );

}

std::string
CustomRealValueMetric::name() const {
	return name_static();
}

std::string
CustomRealValueMetric::name_static() {
	return "CustomRealValueMetric";

}
std::string
CustomRealValueMetric::metric() const {
	return get_custom_type().empty() ? "custom_real_valued_metric" : get_custom_type();
}

void
CustomRealValueMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  )
{
	SimpleMetric::parse_base_tag( tag );
	set_value( tag->getOption< core::Real >( "value", 0.0 ) );
}

void
CustomRealValueMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("value", xsct_real, "The custom floating-point value to cache in the pose.", "");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for caching an arbitrary, user-defined floating-point value in a pose.", attlist);
}

/// @brief Provide the citation.
void
CustomRealValueMetric::provide_citation_info(
	basic::citation_manager::CitationCollectionList & citations
) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"CustomRealValueMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Vikram K. Mulligan",
		"Systems Biology Group, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org"
		)
	);
}

core::Real
CustomRealValueMetric::calculate(const core::pose::Pose & ) const {
	return value_;
}

void
CustomRealValueMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	CustomRealValueMetric::provide_xml_schema( xsd );
}

std::string
CustomRealValueMetricCreator::keyname() const {
	return CustomRealValueMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
CustomRealValueMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< CustomRealValueMetric >();
}

} //core
} //simple_metrics
} //metrics


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::metrics::CustomRealValueMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	arc( CEREAL_NVP( value_ ) );

}

template< class Archive >
void
core::simple_metrics::metrics::CustomRealValueMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );
	arc( value_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::metrics::CustomRealValueMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::metrics::CustomRealValueMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_metrics_CustomRealValueMetric )
#endif // SERIALIZATION





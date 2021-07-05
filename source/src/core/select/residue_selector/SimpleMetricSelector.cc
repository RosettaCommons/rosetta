// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SimpleMetricSelector.hh
/// @brief  Allows selecting residues based on the result of a PerResidueRealSimpleMetric
/// @author Brian Coventry (bcov@uw.edu)

// Unit headers
#include <core/select/residue_selector/SimpleMetricSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.select.residue_selector.SimpleMetricSelector" );

namespace core {
namespace select {
namespace residue_selector {

SimpleMetricSelector::SimpleMetricSelector(
	core::simple_metrics::PerResidueRealMetricCOP metric /*= nullptr*/,
	Real lower_bound /*= utility::get_undefined_real()*/,
	Real upper_bound /*= utility::get_undefined_real()*/,
	bool outside_bounds /*= false*/
) : core::select::residue_selector::ResidueSelector(),
	lower_bound_( lower_bound ),
	upper_bound_( upper_bound ),
	outside_bounds_( outside_bounds )
{
	set_metric( metric );
}

SimpleMetricSelector::~SimpleMetricSelector() = default;

SimpleMetricSelector::SimpleMetricSelector( SimpleMetricSelector const & ot ) :
	core::select::residue_selector::ResidueSelector( ot )
{
	*this = ot;
}

SimpleMetricSelector &
SimpleMetricSelector::operator=( SimpleMetricSelector const & ot ) {
	core::select::residue_selector::ResidueSelector::operator=( ot );
	if ( ot.metric_ ) {
		metric_ = std::dynamic_pointer_cast< core::simple_metrics::PerResidueRealMetric const >( ot.metric_->clone() );
	}
	lower_bound_ = ot.lower_bound_;
	upper_bound_ = ot.upper_bound_;
	outside_bounds_ = ot.outside_bounds_;
	return *this;
}

ResidueSelectorOP
SimpleMetricSelector::clone() const {
	return utility::pointer::make_shared< SimpleMetricSelector >( *this );
}


ResidueSubset
SimpleMetricSelector::apply( core::pose::Pose const & pose ) const {

	if ( ! metric_ ) {
		utility_exit_with_message("SimpleMetricSelector: You must set a SimpleMetric!");
	}

	bool lb_defined = !utility::isnan( lower_bound_ );
	bool ub_defined = !utility::isnan( upper_bound_ );

	if ( ! ( lb_defined || ub_defined ) ) {
		utility_exit_with_message("SimpleMetricSelector: You must define a lower_bound or an upper_bound!");
	}

	if ( outside_bounds_ && ! ( lb_defined && ub_defined ) ) {
		utility_exit_with_message("SimpleMetricSelector: With outside_bounds, you must define both a lower_bound and an upper_bound!");
	}

	ResidueSubset subset( pose.size(), false );

	std::map< Size, Real > metric_result = metric_->calculate( pose );

	for ( auto const & pair : metric_result ) {
		Size seqpos = pair.first;
		Real value = pair.second;

		bool within_lb = lb_defined ? value >= lower_bound_ : true;
		bool within_ub = ub_defined ? value <= upper_bound_ : true;

		if ( outside_bounds_ ) {
			subset[ seqpos ] = ! within_lb || ! within_ub;
		} else {
			subset[ seqpos ] = within_lb && within_ub;
		}
	}


	return subset;
}


void
SimpleMetricSelector::set_metric( core::simple_metrics::PerResidueRealMetricCOP const & metric ) {
	if ( ! metric ) {
		metric_ = nullptr;
	} else {
		metric_ = std::dynamic_pointer_cast< core::simple_metrics::PerResidueRealMetric const >( metric->clone() );
	}
}

bool has_option( utility::tag::TagCOP const & tag, std::string const & key ) {
	std::string check = "hasOptionIsBrokenForEmptyStrings";
	return tag->getOption<std::string>( key, check ) != check;
}

void
SimpleMetricSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->hasOption( "metric" ) ) {
		core::simple_metrics::SimpleMetricCOP metric = core::simple_metrics::get_metric_from_datamap_and_subtags( tag, datamap, "metric" );
		core::simple_metrics::PerResidueRealMetricCOP casted_metric =
			std::dynamic_pointer_cast< core::simple_metrics::PerResidueRealMetric const >( metric );
		if ( ! casted_metric ) {
			utility_exit_with_message("SimpleMetricSelector: You must pass a PerResidueRealMetric to this!");
		}
		set_metric( casted_metric );
	}
	// Leaving lower_bound blank is valid
	if ( has_option( tag, "lower_bound" ) && !tag->getOption<std::string>( "lower_bound" ).empty() ) {
		lower_bound_ = tag->getOption<Real>( "lower_bound" );
	}
	// Leaving upper_bound blank is valid
	if ( has_option( tag, "upper_bound" ) && !tag->getOption<std::string>( "upper_bound" ).empty() ) {
		upper_bound_ = tag->getOption<Real>( "upper_bound" );
	}
	if ( tag->hasOption( "outside_bounds" ) ) {
		outside_bounds_ = tag->getOption<bool>( "outside_bounds" );
	}
}


std::string SimpleMetricSelector::get_name() const {
	return SimpleMetricSelector::class_name();
}

std::string SimpleMetricSelector::class_name() {
	return "SimpleMetricSelector";
}

void
SimpleMetricSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute(
		"metric", xs_string ,
		"A PerResidueRealMetric that will be used to select residues." )
		+ XMLSchemaAttribute::attribute_w_default(
		"lower_bound", xs_string,
		"The lowest value that will be considered \"inside the bounds\". Leaving this blank or not including it will have no lower_bound.",
		"" )
		+ XMLSchemaAttribute::attribute_w_default(
		"upper_bound", xs_string,
		"The highest value that will be considered \"inside the bounds\". Leaving this blank or not including it will have no upper_bound.",
		"" )
		+ XMLSchemaAttribute::attribute_w_default(
		"outside_bounds", xsct_rosetta_bool,
		"Use this if you want to select residues below lower_bound and above upper_bound.",
		"false" );

	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(),
		"Allows selecting residues based on a PerResidueRealMetric.", attlist );
}


ResidueSelectorOP
SimpleMetricSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared<SimpleMetricSelector>();
}

std::string
SimpleMetricSelectorCreator::keyname() const {
	return SimpleMetricSelector::class_name();
}

void
SimpleMetricSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SimpleMetricSelector::provide_xml_schema( xsd );
}


} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::SimpleMetricSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( metric_ ) );
	arc( CEREAL_NVP( lower_bound_ ) );
	arc( CEREAL_NVP( upper_bound_ ) );
	arc( CEREAL_NVP( outside_bounds_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::SimpleMetricSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	std::shared_ptr< core::simple_metrics::PerResidueRealMetric > local_metric;
	arc( local_metric );
	metric_ = local_metric; // copy the non-const pointer(s) into the const pointer(s)
	arc( lower_bound_ );
	arc( upper_bound_ );
	arc( outside_bounds_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::SimpleMetricSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::SimpleMetricSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_SimpleMetricSelector )
#endif // SERIALIZATION

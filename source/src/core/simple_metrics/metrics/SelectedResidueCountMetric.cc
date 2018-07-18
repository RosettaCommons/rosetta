// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SelectedResidueCountMetric.cc
/// @brief A SimpleMetric that counts the number of residues in a residue selection.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <core/simple_metrics/metrics/SelectedResidueCountMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.SelectedResidueCountMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SelectedResidueCountMetric::SelectedResidueCountMetric():
	core::simple_metrics::RealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SelectedResidueCountMetric::~SelectedResidueCountMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SelectedResidueCountMetric::SelectedResidueCountMetric( SelectedResidueCountMetric const & src ):
	core::simple_metrics::RealMetric( src ),
	residue_selector_(nullptr)
{

}


core::simple_metrics::SimpleMetricOP
SelectedResidueCountMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new SelectedResidueCountMetric( *this ) );
}

std::string
SelectedResidueCountMetric::name() const {
	return name_static();
}

std::string
SelectedResidueCountMetric::name_static() {
	return "SelectedResidueCountMetric";

}
std::string
SelectedResidueCountMetric::metric() const {
	return "selection_count";
}

void
SelectedResidueCountMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &datamap  )
{
	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector( core::select::residue_selector::parse_residue_selector(tag, datamap, "residue_selector" ) );
	}
}

void
SelectedResidueCountMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attributes_for_parse_residue_selector( attlist, "residue_selector", "A residue selector.  The number of residues selected by this selector will be returned as the count.  If not provided, the number of residues in the pose will be returned." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for counting the number of residues in a pose or residue selection, and adding the count to the resulting score file.", attlist);
}

core::Real
SelectedResidueCountMetric::calculate(const core::pose::Pose &pose ) const {
	if ( residue_selector_ == nullptr ) return static_cast< core::Real >(pose.total_residue());
	utility::vector1< bool > selection( residue_selector_->apply(pose) );
	core::Size count(0);
	for ( core::Size i(1), imax(selection.size()); i<=imax; ++i ) {
		if ( selection[i] ) ++count;
	}
	return static_cast< core::Real >( count );
}

/// @brief Set the residue selector.
/// @details  Copies the input pointer; doesn't clone the object.
void
SelectedResidueCountMetric::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	runtime_assert_string_msg( selector_in != nullptr, "Error in SelectedResidueCountMetric::set_residue_selector(): The pointer passed to this function was null.");
	residue_selector_ = selector_in;
}

/// @brief Remove the residue selector.
/// @details In the absence of a residue selector, this metric returns the number of residues in a pose.
void
SelectedResidueCountMetric::remove_residue_selector() {
	residue_selector_ = nullptr;
}

void
SelectedResidueCountMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SelectedResidueCountMetric::provide_xml_schema( xsd );
}

std::string
SelectedResidueCountMetricCreator::keyname() const {
	return SelectedResidueCountMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SelectedResidueCountMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new SelectedResidueCountMetric );

}

} //protocols
} //analysis
} //simple_metrics







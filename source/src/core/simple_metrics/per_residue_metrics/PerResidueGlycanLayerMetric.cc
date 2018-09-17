// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueGlycanLayerMetric.cc
/// @brief A metric that outputs the layer of the glycan tree as measured by the residue distance to the root.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/PerResidueGlycanLayerMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.PerResidueGlycanLayerMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select;
using namespace core::select::residue_selector;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PerResidueGlycanLayerMetric::PerResidueGlycanLayerMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueGlycanLayerMetric::~PerResidueGlycanLayerMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PerResidueGlycanLayerMetric::PerResidueGlycanLayerMetric( PerResidueGlycanLayerMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
PerResidueGlycanLayerMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new PerResidueGlycanLayerMetric( *this ) );

}

std::string
PerResidueGlycanLayerMetric::name() const {
	return name_static();
}

std::string
PerResidueGlycanLayerMetric::name_static() {
	return "PerResidueGlycanLayerMetric";

}
std::string
PerResidueGlycanLayerMetric::metric() const {
	return "glycan_layer";
}

void
PerResidueGlycanLayerMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );

}

void
PerResidueGlycanLayerMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	//attributes_for_parse_residue_selector( attlist, "residue_selector",
	// "Selector specifying residues." );

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for getting the layer of glycan residues.  The layer is defined as the residue distance to the start of the glycan tree", attlist);
}

std::map< core::Size, core::Real >
PerResidueGlycanLayerMetric::calculate(const pose::Pose & pose) const {
	utility::vector1< bool > subset = get_selector()->apply(pose);
	std::map< core::Size, core::Real > result;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! subset[i] ) continue;
		if ( pose.residue_type( i ).is_carbohydrate() ) {
			core::Size layer = pose.glycan_tree_set()->get_distance_to_start( i );
			result[i] = layer;
		}
	}
	return result;

}

void
PerResidueGlycanLayerMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueGlycanLayerMetric::provide_xml_schema( xsd );
}

std::string
PerResidueGlycanLayerMetricCreator::keyname() const {
	return PerResidueGlycanLayerMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueGlycanLayerMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new PerResidueGlycanLayerMetric );

}

} //core
} //simple_metrics
} //per_residue_metrics







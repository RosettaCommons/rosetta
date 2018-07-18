// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueSasaMetric.cc
/// @brief A per-residue metric that will calculate SASA for each residue given in a selector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/PerResidueSasaMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/scoring/sasa/SasaCalc.hh>
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


static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.PerResidueSasaMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select::residue_selector;
using namespace core::scoring::sasa;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PerResidueSasaMetric::PerResidueSasaMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueSasaMetric::~PerResidueSasaMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PerResidueSasaMetric::PerResidueSasaMetric( PerResidueSasaMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
PerResidueSasaMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new PerResidueSasaMetric( *this ) );

}

std::string
PerResidueSasaMetric::name() const {
	return name_static();
}

std::string
PerResidueSasaMetric::name_static() {
	return "PerResidueSasaMetric";

}
std::string
PerResidueSasaMetric::metric() const {
	return "res_sasa";
}


void
PerResidueSasaMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );

}

void
PerResidueSasaMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for measuring per-residue SASA and adding it to the resulting score file..", attlist);
}

std::map< core::Size, core::Real >
PerResidueSasaMetric::calculate(const pose::Pose & pose) const {
	SasaCalc calc = SasaCalc();
	calc.calculate( pose );

	std::map< core::Size, core::Real > result;
	utility::vector1< core::Real > rsd_sasa = calc.get_residue_sasa();
	for ( core::Size res : core::select::get_residues_from_subset( get_selector()->apply( pose )) ) {
		result[res] = rsd_sasa[res];
	}
	return result;
}

void
PerResidueSasaMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueSasaMetric::provide_xml_schema( xsd );
}

std::string
PerResidueSasaMetricCreator::keyname() const {
	return PerResidueSasaMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueSasaMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new PerResidueSasaMetric );

}

} //core
} //simple_metrics
} //per_residue_metrics







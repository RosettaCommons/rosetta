// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/metrics/SasaMetric.cc
/// @brief A Metric to cacluate overall sasa of a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/simple_metrics/metrics/SasaMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueSasaMetric.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "core.simple_metrics.metrics.SasaMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

using namespace core::scoring::sasa;
using namespace core::select::residue_selector;

SasaMetric::SasaMetric():
	RealMetric()
{

}

SasaMetric::SasaMetric( select::residue_selector::ResidueSelectorCOP selector):
	RealMetric()
{
	set_residue_selector( selector );
}

SasaMetric::~SasaMetric(){}

SasaMetric::SasaMetric( SasaMetric const & src):
	RealMetric( src )
{
	selector_ = src.selector_;
}


std::string
SasaMetric::name() const {
	return name_static();
}

std::string
SasaMetric::name_static() {
	return "SasaMetric";

}
std::string
SasaMetric::metric() const {
	return "sasa";
}

void
SasaMetric::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

core::Real
SasaMetric::calculate(const pose::Pose &pose) const {
	using namespace core::simple_metrics::per_residue_metrics;

	PerResidueSasaMetric core_metric = PerResidueSasaMetric();

	if ( selector_ ) {
		core_metric.set_residue_selector( selector_);

	}
	std::map< core::Size, core::Real > res_sasa = core_metric.calculate( pose );
	core::Real total_sasa = 0;
	for ( auto res_sasa_pair : res_sasa ) {
		total_sasa+= res_sasa_pair.second;
	}
	return total_sasa;
}

void
SasaMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap  )
{

	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(select::residue_selector::parse_residue_selector( tag, datamap ));
	}

}

SimpleMetricOP
SasaMetric::clone() const {
	return SimpleMetricOP(new SasaMetric( *this ) );

}

void
SasaMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	std::string description = "If a residue selector is present, we calculate the total sasa of these residues." ;

	attributes_for_parse_residue_selector_default_option_name(attlist, description );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for measuring SASA and adding it to the resulting score file.", attlist);
}

void
SasaMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SasaMetric::provide_xml_schema( xsd );
}

std::string
SasaMetricCreator::keyname() const {
	return SasaMetric::name_static();
}

SimpleMetricOP
SasaMetricCreator::create_simple_metric() const {
	return SimpleMetricOP( new SasaMetric );

}

} //core
} //simple_metrics
} //metrics







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
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added support for polar or hydrophobic SASA.

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

SasaMetric::SasaMetric(
	select::residue_selector::ResidueSelectorCOP selector,
	core::scoring::sasa::SasaMethodHPMode const mode /*= core::scoring::sasa::SasaMethodHPMode::ALL_SASA*/
):
	RealMetric()
{
	set_residue_selector( selector );
	set_sasa_metric_mode( mode );
}

SasaMetric::~SasaMetric(){}


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

/// @brief Set the behaviour of this metric (count all SASA, count polar SASA, count hydrophobic SASA, etc.).
/// @details Default is to compute all SASA.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
SasaMetric::set_sasa_metric_mode(
	core::scoring::sasa::SasaMethodHPMode const mode_in
) {
	runtime_assert_string_msg( mode_in < core::scoring::sasa::SasaMethodHPMode::END_OF_LIST && static_cast< core::Size >(mode_in) > 0, "Error in SasaMetric::set_sasa_metric_mode(): Unrecognized mode provided!" );
	sasa_metric_mode_ = mode_in;
}

/// @brief Set the behaviour of this metric (count all SASA, count polar SASA, count hydrophobic SASA, etc.).
/// @details Default is to compute all SASA.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
SasaMetric::set_sasa_metric_mode(
	std::string const & mode_in
) {
	core::scoring::sasa::SasaMethodHPMode const mode( core::scoring::sasa::SasaMethod::sasa_metric_mode_from_name( mode_in ) );
	runtime_assert_string_msg( mode != core::scoring::sasa::SasaMethodHPMode::INVALID_MODE, "Error in SasaMetric::set_sasa_metric_mode(): Could not parse \"" + mode_in + "\" as a valid mode." );
	set_sasa_metric_mode( mode );
}

void
SasaMetric::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

core::Real
SasaMetric::calculate(const pose::Pose &pose) const {
	using namespace core::simple_metrics::per_residue_metrics;

	PerResidueSasaMetric core_metric( sasa_metric_mode_ );

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

	set_sasa_metric_mode( tag->getOption<std::string>("sasa_metric_mode", "all_sasa") );
}

SimpleMetricOP
SasaMetric::clone() const {
	return utility::pointer::make_shared< SasaMetric >( *this );

}

void
SasaMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	std::string description = "If a residue selector is present, we calculate the total sasa of these residues." ;

	attributes_for_parse_residue_selector_default_option_name(attlist, description );

	attlist + XMLSchemaAttribute::attribute_w_default( "sasa_metric_mode", xs_string, "Sets the behaviour of the calculator (the subset of the SASA that is counted).  Options include: " + core::scoring::sasa::SasaMethod::list_sasa_method_hp_modes() + ".", "all_sasa" );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"A metric for measuring SASA and adding it to the resulting score file.  Modified 19 Aug. 2019 by Vikram K. Mulligan to add options for polar or hydrophobic SASA.", attlist);
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
	return utility::pointer::make_shared< SasaMetric >();

}

} //core
} //simple_metrics
} //metrics







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
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added support for polar or hydrophobic

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

/// @brief Mode constructor.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
PerResidueSasaMetric::PerResidueSasaMetric(
	core::scoring::sasa::SasaMethodHPMode const mode
) {
	set_mode( mode );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueSasaMetric::~PerResidueSasaMetric(){}

core::simple_metrics::SimpleMetricOP
PerResidueSasaMetric::clone() const {
	return utility::pointer::make_shared< PerResidueSasaMetric >( *this );

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

	set_mode( tag->getOption<std::string>("mode", "all_sasa") );
}

void
PerResidueSasaMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;
	using namespace core::scoring::sasa;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "mode", xs_string, "Sets the behaviour of the calculator (the subset of the SASA that is counted).  Options include: " + SasaMethod::list_sasa_method_hp_modes() + ".", "all_sasa" );

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"A metric for measuring per-residue SASA and adding it to the resulting score file.  Modified 19 Aug. 2019 by Vikram K. Mulligan to add options for polar or hydrophobic SASA.", attlist);
}

std::map< core::Size, core::Real >
PerResidueSasaMetric::calculate(const pose::Pose & pose) const {
	SasaCalc calc;
	calc.set_sasa_method_hp_mode( mode_ );
	calc.calculate( pose );

	std::map< core::Size, core::Real > result;
	utility::vector1< core::Real > rsd_sasa = calc.get_residue_sasa();
	for ( core::Size res : core::select::get_residues_from_subset( get_selector()->apply( pose )) ) {
		result[res] = rsd_sasa[res];
	}
	return result;
}

/// @brief Set the SASA mode (all SASA, polar only, hydrophobic only, etc.).
/// @note The enum class is defined in core/simple_metrics/metrics/SasaMetric.hh.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PerResidueSasaMetric::set_mode(
	core::scoring::sasa::SasaMethodHPMode const mode_in
) {
	using namespace core::scoring::sasa;
	runtime_assert_string_msg( mode_in < SasaMethodHPMode::END_OF_LIST && static_cast< core::Size >(mode_in) > 0, "Error in PerResidueSasaMetric::set_mode(): Unrecognized mode provided!" );
	mode_ = mode_in;
}

/// @brief Set the SASA mode (all SASA, polar only, hydrophobic only, etc.).
/// @note The enum class is defined in core/simple_metrics/metrics/SasaMetric.hh.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PerResidueSasaMetric::set_mode(
	std::string const & mode_in
) {
	using namespace core::scoring::sasa;
	SasaMethodHPMode const mode( SasaMethod::sasa_metric_mode_from_name( mode_in ) );
	runtime_assert_string_msg( mode != SasaMethodHPMode::INVALID_MODE, "Error in PerResidueSasaMetric::set_mode(): Could not parse \"" + mode_in + "\" as a valid mode." );
	set_mode( mode );
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
	return utility::pointer::make_shared< PerResidueSasaMetric >();

}

} //core
} //simple_metrics
} //per_residue_metrics







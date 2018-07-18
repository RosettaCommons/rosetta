// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SecondaryStructureMetric.cc
/// @brief A SimpleMetric to output the secondary structure of the protein or residues selected by a residue selector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/SecondaryStructureMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.SecondaryStructureMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SecondaryStructureMetric::SecondaryStructureMetric():
	core::simple_metrics::StringMetric()
{}

SecondaryStructureMetric::SecondaryStructureMetric( select::residue_selector::ResidueSelectorCOP selector):
	core::simple_metrics::StringMetric()
{
	set_residue_selector( selector );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SecondaryStructureMetric::~SecondaryStructureMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SecondaryStructureMetric::SecondaryStructureMetric( SecondaryStructureMetric const & src ):
	core::simple_metrics::StringMetric( src ),
	dssp_reduced_(src.dssp_reduced_)
{
	selector_ = src.selector_;
}


core::simple_metrics::SimpleMetricOP
SecondaryStructureMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new SecondaryStructureMetric( *this ) );

}

std::string
SecondaryStructureMetric::name() const {
	return name_static();
}

std::string
SecondaryStructureMetric::name_static() {
	return "SecondaryStructureMetric";

}
std::string
SecondaryStructureMetric::metric() const {
	return "secondary_structure";
}

void
SecondaryStructureMetric::set_residue_selector(select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
SecondaryStructureMetric::set_use_dssp_reduced(bool reduced){
	dssp_reduced_ = reduced;
}

void
SecondaryStructureMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap)
{
	SimpleMetric::parse_base_tag( tag );

	set_use_dssp_reduced( tag->getOption< bool >("dssp_reduced", true));

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(core::select::residue_selector::parse_residue_selector( tag, datamap ));
	}

}

void
SecondaryStructureMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("dssp_reduced", xsct_rosetta_bool, "Use the reduced DSSP alphabet (3 letters)", "true");

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Use a residue selector to output the secondary structure of only residues in selection." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric getting the Secondary Structure string of a pose or subset of residues via residue selector using DSSP", attlist);
}

std::string
SecondaryStructureMetric::calculate(const pose::Pose & pose) const {

	utility::vector1< bool > subset( pose.size(), true);
	if ( selector_ ) {
		subset = selector_->apply( pose );
	}
	std::string sec_struct = "";
	utility::vector1< core::Size > selection = select::get_residues_from_subset( subset );

	core::scoring::dssp::Dssp dssp( pose );

	if ( dssp_reduced_ ) dssp.dssp_reduced();

	for ( core::Size resnum : selection ) {
		sec_struct += utility::to_string( dssp.get_dssp_secstruct( resnum ) );
	}

	return sec_struct;
}

void
SecondaryStructureMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SecondaryStructureMetric::provide_xml_schema( xsd );
}

std::string
SecondaryStructureMetricCreator::keyname() const {
	return SecondaryStructureMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SecondaryStructureMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new SecondaryStructureMetric );

}

} //core
} //simple_metrics
} //metrics







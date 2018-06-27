// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SequenceMetric.cc
/// @brief A SimpleMetric to output the single-letter OR three-letter sequence of a protein or subset of positions/regions using a ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/SequenceMetric.hh>
#include <core/simple_metrics/metrics/SequenceMetricCreator.hh>

// Core headers
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.SequenceMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SequenceMetric::SequenceMetric():
	core::simple_metrics::StringMetric()
{}

SequenceMetric::SequenceMetric( select::residue_selector::ResidueSelectorCOP selector ):
	core::simple_metrics::StringMetric()
{
	set_residue_selector( selector );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SequenceMetric::~SequenceMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SequenceMetric::SequenceMetric( SequenceMetric const & src ):
	core::simple_metrics::StringMetric( src ),
	three_letter_(src.three_letter_)
{
	selector_ = src.selector_;
}

core::simple_metrics::SimpleMetricOP
SequenceMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new SequenceMetric( *this ) );

}

std::string
SequenceMetric::name() const {
	return name_static();
}

std::string
SequenceMetric::name_static() {
	return "SequenceMetric";

}
std::string
SequenceMetric::metric() const {
	return "sequence";
}

void
SequenceMetric::set_residue_selector(select::residue_selector::ResidueSelectorCOP selector ){
	selector_ = selector;
}

void
SequenceMetric::set_use_three_letter_code(bool three_letter){
	three_letter_ = three_letter;
}

void
SequenceMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap)
{
	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(core::select::residue_selector::parse_residue_selector( tag, datamap ));
	}
	set_use_three_letter_code( tag->getOption< bool >("three_letter", false));
}

void
SequenceMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("three_letter", xsct_rosetta_bool, "Set to ouput the sequence as three-letter codes. Useful for modifications and glycans.  Comma-separated", "false");

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Output the sequence of only the selected residues." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for measuring ... and adding it to the resulting score file.", attlist);
}

std::string
SequenceMetric::calculate(const pose::Pose & pose ) const {

	utility::vector1< bool > subset( pose.size(), true);
	if ( selector_ ) {
		subset = selector_->apply( pose );
	}

	std::string sequence = "";
	utility::vector1< core::Size > selection = select::get_residues_from_subset( subset );
	//std::cout << "Selection"  << utility::to_string(selection);

	for ( core::Size resnum : selection ) {
		if ( three_letter_ ) {
			std::string three = pose.residue_type( resnum ).name3();
			if ( resnum != selection.size() ) {
				three = three + ",";
			}
			sequence += three;
		} else {
			//std::cout << pose.residue_type( resnum ).name1() << std::endl;
			sequence += pose.residue_type( resnum ).name1();
		}
	}
	return sequence;
}

void
SequenceMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SequenceMetric::provide_xml_schema( xsd );
}

std::string
SequenceMetricCreator::keyname() const {
	return SequenceMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SequenceMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new SequenceMetric );

}

} //core
} //simple_metrics
} //metrics







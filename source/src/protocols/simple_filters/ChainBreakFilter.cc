// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ChainBreakFilter.cc
/// @brief
/// @author Christoffer Norn


//Unit Headers
#include <protocols/simple_filters/ChainBreakFilter.hh>
#include <protocols/simple_filters/ChainBreakFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
//#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.ChainBreak" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP ChainBreakFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ChainBreak ); }

// XRW TEMP std::string
// XRW TEMP ChainBreakFilterCreator::keyname() const { return "ChainBreak"; }

//default ctor
ChainBreak::ChainBreak() :
	protocols::filters::Filter( "ChainBreak" ),
	threshold_( 1 ),
	chain_num_( 1 ),
	tolerance_( 0.13 )
{
}

ChainBreak::~ChainBreak() = default;

void
ChainBreak::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	threshold( tag->getOption< core::Size >( "threshold", 1 ) );
	chain_num( tag->getOption< core::Size >( "chain_num", 1 ) );
	tolerance( tag->getOption< core::Real >( "tolerance", 0.13 ) );

	TR << "ChainBreak with options: threshold: " << threshold() << " chain_num: "<< chain_num() << std::endl;
}

bool
ChainBreak::apply( core::pose::Pose const & pose ) const {
	core::Size val = compute( pose );
	TR << "ChainBreak filter identified " << val << " chainbreaks." << std::endl;
	return( compute( pose ) <= threshold() );
}

void
ChainBreak::report( std::ostream & o, core::pose::Pose const & pose ) const {
	bool const val = ( compute( pose ) >= threshold() );
	o << "ChainBreak returns " << val << std::endl;
}

core::Real
ChainBreak::report_sm( core::pose::Pose const & pose ) const {
	return( core::Real (compute( pose )) );
}

core::Size
ChainBreak::compute(
	core::pose::Pose const & pose
) const {
	core::Size chainBreaks = 0;
	core::Real max_bl = 1.33 + tolerance();
	core::Real min_bl = 1.33 - tolerance();
	TR << "Will check peptide bond lengths between " << pose.conformation().chain_begin( chain_num() ) << " to " << pose.conformation().chain_end( chain_num() )-1 << std::endl;
	TR<< "bond length tolerance value is:" << tolerance() << std::endl;
	for ( core::Size resj = pose.conformation().chain_begin( chain_num() ); resj <= pose.conformation().chain_end( chain_num() ) - 1; ++resj ) {
		core::Real const distance = pose.residue( resj + 1 ).xyz( "N" ).distance(pose.residue( resj ).xyz( "C" ));
		if ( distance > max_bl || distance < min_bl ) {
			TR<<"distance is: "<<distance<<" max_bl"<<max_bl<<std::endl;
			TR << "The distance from " << resj << " to " << resj+1 << " is " << distance << std::endl;
			++chainBreaks;
		}
	}
	return chainBreaks;
}

std::string ChainBreak::name() const {
	return class_name();
}

std::string ChainBreak::class_name() {
	return "ChainBreak";
}

void ChainBreak::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"threshold", xsct_non_negative_integer,
		"Number of chainbreaks allowed",
		"1")
		+ XMLSchemaAttribute::attribute_w_default(
		"chain_num", xsct_non_negative_integer,
		"which chain should we check for",
		"1")
		+ XMLSchemaAttribute::attribute_w_default(
		"tolerance", xsct_real,
		"the allowed angstrom deviation from the mean optimal bond length",
		"0.13");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Measures the number of chainBreaks in the pose",
		attlist );
}


std::string ChainBreakFilterCreator::keyname() const {
	return ChainBreak::class_name();
}

protocols::filters::FilterOP
ChainBreakFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ChainBreak );
}

void ChainBreakFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ChainBreak::provide_xml_schema( xsd );
}


}
}

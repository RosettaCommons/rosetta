// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)

#include <protocols/filters/ContingentFilter.hh>
#include <protocols/filters/ContingentFilterCreator.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace filters {

static basic::Tracer TR( "protocols.filters.ContingentFilter" );

/// @brief default ctor
ContingentFilter::ContingentFilter() :
	parent( "ContingentFilter" ),
	value_( false )
{}

bool
ContingentFilter::apply(core::pose::Pose const & ) const
{
	return( get_value() );
}

core::Real
ContingentFilter::report_sm( core::pose::Pose const & ) const
{
	return( get_value() );
}

void
ContingentFilter::report( std::ostream & out, core::pose::Pose const & ) const
{
	out<<"ContingentFilter returns "<<get_value()<<std::endl;
}

void
ContingentFilter::set_value( bool const value ){
	value_ = value;
}

bool
ContingentFilter::get_value() const{
	return( value_ );
}

void
ContingentFilter::parse_my_tag( utility::tag::TagCOP const,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &)
{
	TR.Info << "ContingentFilter"<<std::endl;
}

// XRW TEMP FilterOP
// XRW TEMP ContingentFilterCreator::create_filter() const { return FilterOP( new ContingentFilter ); }

// XRW TEMP std::string
// XRW TEMP ContingentFilterCreator::keyname() const { return "ContingentFilter"; }

std::string ContingentFilter::name() const {
	return class_name();
}

std::string ContingentFilter::class_name() {
	return "ContingentFilter";
}

void ContingentFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"A special filter that allows movers to set its value (pass/fail). "
		"This value can then be used in the protocol together with IfMover to "
		"control the flow of execution depending on the success of the mover.",
		attlist );
}

std::string ContingentFilterCreator::keyname() const {
	return ContingentFilter::class_name();
}

protocols::filters::FilterOP
ContingentFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ContingentFilter );
}

void ContingentFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ContingentFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ExpiryFilter.cc
/// @brief
/// @author Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/ExpiryFilter.hh>
#include <protocols/simple_filters/ExpiryFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <ctime>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.simple_filters.ExpiryFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP ExpiryFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ExpiryFilter ); }

// XRW TEMP std::string
// XRW TEMP ExpiryFilterCreator::keyname() const { return "Expiry"; }

//default ctor
ExpiryFilter::ExpiryFilter() :
	protocols::filters::Filter( "Expiry" ),
	seconds_( 0 )
{
	start_time_ = time( nullptr );
}

ExpiryFilter::~ExpiryFilter() = default;

void
ExpiryFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	seconds( tag->getOption< core::Size >( "seconds" ) );
}

bool
ExpiryFilter::apply( core::pose::Pose const & pose ) const {
	return( compute( pose ) <= seconds() );
}

void
ExpiryFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out<<"Time from start "<< compute( pose )<<" seconds"<<std::endl;
}

core::Real
ExpiryFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
ExpiryFilter::compute(
	core::pose::Pose const & /*pose*/
) const {
	return( time( nullptr ) - start_time() );
}

core::Size
ExpiryFilter::seconds() const{
	return seconds_;
}

void
ExpiryFilter::seconds( core::Size const s ){
	seconds_ = s;
}

core::Size
ExpiryFilter::start_time() const{
	return start_time_;
}

std::string ExpiryFilter::name() const {
	return class_name();
}

std::string ExpiryFilter::class_name() {
	return "Expiry";
}

void ExpiryFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("seconds", xsct_non_negative_integer, "how many seconds until this triggers failure?");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Has a predetermined number of seconds elapsed since the start of the trajectory? If so, return false (to stop the trajectory), else return true.", attlist );
}

std::string ExpiryFilterCreator::keyname() const {
	return ExpiryFilter::class_name();
}

protocols::filters::FilterOP
ExpiryFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ExpiryFilter );
}

void ExpiryFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExpiryFilter::provide_xml_schema( xsd );
}


}
}

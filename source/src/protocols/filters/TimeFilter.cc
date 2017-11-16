// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/TimeFilter.cc
/// @author Sarel Fleishman (sarelf@uw.edu)

#include <protocols/filters/TimeFilter.hh>
#include <protocols/filters/TimeFilterCreator.hh>

// Project Headers
#include <core/types.hh>

//parsing
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace filters {

static basic::Tracer TR( "protocols.filters.TimeFilter" );

/// @brief default ctor
TimeFilter::TimeFilter() :
	Filter( "Time" ),
	tic_( 0 ),
	toc_( 0 )
{}

TimeFilter::~TimeFilter() = default;

core::Real
TimeFilter::compute( core::pose::Pose const & ) const
{
	if ( tic_ == 0 ) {
		tic_ = time( nullptr );
	} else {
		toc_ = time( nullptr );
		return toc_ - tic_;
	}
	return 0;
}

core::Real
TimeFilter::report_sm( core::pose::Pose const & ) const
{
	core::Real const elapsed_time( toc_ - tic_ );
	return( elapsed_time );
}

void
TimeFilter::report( std::ostream &out, core::pose::Pose const & ) const
{
	core::Real const elapsed_time( toc_ - tic_ );
	out<<elapsed_time<<" seconds"<<std::endl;
}

void TimeFilter::parse_my_tag( utility::tag::TagCOP const,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &)
{
}

bool
TimeFilter::apply( core::pose::Pose const & p ) const{
	compute( p );
	return true;
}

// XRW TEMP FilterOP
// XRW TEMP TimeFilterCreator::create_filter() const { return FilterOP( new TimeFilter ); }

// XRW TEMP std::string
// XRW TEMP TimeFilterCreator::keyname() const { return "Time"; }

std::string TimeFilter::name() const {
	return class_name();
}

std::string TimeFilter::class_name() {
	return "Time";
}

void TimeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Report how long a sequence of movers/filters takes. Within the "
		"protocol, you need to call time at least twice, once, when you "
		"want to start the timer, and then, when you want to stop. "
		"The reported time is that between the first and last calls.",
		attlist );
}

std::string TimeFilterCreator::keyname() const {
	return TimeFilter::class_name();
}

protocols::filters::FilterOP
TimeFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new TimeFilter );
}

void TimeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TimeFilter::provide_xml_schema( xsd );
}


} // filters
} // protocols

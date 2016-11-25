// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/NcontactsFilter.cc
/// @brief filter structures by sheet topology
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/NcontactsFilter.hh>
#include <protocols/fldsgn/filters/NcontactsFilterCreator.hh>

// Package Headers
#include <protocols/fldsgn/NcontactsCalculator.hh>

// Project Headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

//// C++ headers
#include <map>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.fldsgn.filters.NcontactsFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

/// @brief default constructor
NcontactsFilter::NcontactsFilter():
	Filter( "Ncontacts" ),
	report_type_( "sidechain_heavy_apolar_atm" ),
	filter_value_( 0.0 )
{}

/// @brief default constructor
NcontactsFilter::NcontactsFilter(
	String const & report_type,
	Real const filter_value ):
	Filter( "Ncontacts" ),
	report_type_( report_type ),
	filter_value_( filter_value )
{}

/// @brief copy constructor
NcontactsFilter::NcontactsFilter( NcontactsFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	report_type_( rval.report_type_ ),
	filter_value_( rval.filter_value_ )
{}

/// @brief destructor
NcontactsFilter::~NcontactsFilter(){}

/// @brief
NcontactsFilter::Real
NcontactsFilter::report_sm( Pose const & pose ) const
{
	return compute( pose );
}

/// @brief
void
NcontactsFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "Ncontacts: " << compute( pose ) << "\n";
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool
NcontactsFilter::apply( Pose const & pose ) const
{
	Real score = compute( pose );

	if ( filter_value_ < score ) {
		return true;
	} else {
		return false;
	}
} // apply

/// @brief comute ncontacts
NcontactsFilter::Real
NcontactsFilter::compute( Pose const & pose ) const
{

	basic::MetricValue< Real > ncontact;

	NcontactsCalculator calc_ncon;
	calc_ncon.get( report_type_, ncontact, pose );

	return ncontact.value();
} // compute

/// @brief parse xml
void
NcontactsFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	report_type_ = tag->getOption< String >( "type", "sidechain_heavy_apolar_atm" );
	filter_value_ = tag->getOption< Real >( "value", 0.0 );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP NcontactsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new NcontactsFilter ); }

// XRW TEMP std::string
// XRW TEMP NcontactsFilterCreator::keyname() const { return "Ncontacts"; }

std::string NcontactsFilter::name() const {
	return class_name();
}

std::string NcontactsFilter::class_name() {
	return "Ncontacts";
}

void NcontactsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction type;
	type.name( "ncontacts_filter_type" );
	type.base_type( xs_string );
	type.add_restriction( xsr_enumeration, "sidechain_heavy_apolar_atm" );
	type.add_restriction( xsr_enumeration, "all_heavy_atm" );
	type.add_restriction( xsr_enumeration, "sidechain_heavy_atm_apolar_aa" );
	type.add_restriction( xsr_enumeration, "ss_entropy" );
	xsd.add_top_level_element( type );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "type", "ncontacts_filter_type", "XRW TO DO", "sidechain_heavy_apolar_atm" )
		+ XMLSchemaAttribute::attribute_w_default( "value", xsct_real, "XRW TO DO", "0.0" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string NcontactsFilterCreator::keyname() const {
	return NcontactsFilter::class_name();
}

protocols::filters::FilterOP
NcontactsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new NcontactsFilter );
}

void NcontactsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NcontactsFilter::provide_xml_schema( xsd );
}



} // filters
} // fldsgn
} // protocols

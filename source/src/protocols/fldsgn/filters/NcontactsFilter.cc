// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/NcontactsFilter.cc
/// @brief filter structures by sheet topology
/// @detailed
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
// AUTO-REMOVED #include <cmath>
#include <map>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


static basic::Tracer tr("protocols.fldsgn.filters.NcontactsFilter");

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

	if( filter_value_ < score ){
		return true;
	}else{
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

protocols::filters::FilterOP
NcontactsFilterCreator::create_filter() const { return new NcontactsFilter; }

std::string
NcontactsFilterCreator::keyname() const { return "Ncontacts"; }



} // filters
} // fldsgn
} // protocols

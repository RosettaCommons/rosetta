// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/RangeFilter.cc
/// @brief
/// @details
/// @author Javier Castellanos (javiercv@uw.edu)

// Unit Headers
#include <protocols/simple_filters/RangeFilter.hh>
#include <protocols/simple_filters/RangeFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

//// C++ headers
static THREAD_LOCAL basic::Tracer tr( "protocols.filters.RangeFilter" );

namespace protocols {
namespace simple_filters {

// @brief default constructor
RangeFilter::RangeFilter() : protocols::filters::Filter()
{}

RangeFilter::RangeFilter(Real lower_bound, Real upper_bound, FilterOP const & filter ) : protocols::filters::Filter(),
	filter_(filter),
	lower_bound_( lower_bound ),
	upper_bound_( upper_bound )
{}

// @brief copy constructor
RangeFilter::RangeFilter( RangeFilter const & rval ) : protocols::filters::Filter(),
	filter_( rval.filter_),
	lower_bound_( rval.lower_bound_ ),
	upper_bound_( rval.upper_bound_ )
{}


/// @brief
void
RangeFilter::report( std::ostream & out, Pose const & pose ) const
{
	Real value = filter_->apply( pose );
	out << value << " in range " << lower_bound_ << " - " << upper_bound_ << std::endl;
}


// @brief returns true if the given pose passes the filter, false otherwise.
bool RangeFilter::apply( Pose const & pose ) const
{
	Real value = filter_->report_sm( pose );
	if ( value > lower_bound_ && value < upper_bound_ ) {
		tr << "Successfully filtered: " << value << " in range " << lower_bound_ << " - " << upper_bound_ << std::endl;
		return true;
	} else {
		tr << "Filter failed: value = " << value << " range = "<< lower_bound_ << " - " << upper_bound_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
RangeFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &filters,
	Movers_map const &,
	Pose const & )
{
	std::string const filter_name( tag->getOption< std::string >( "filter") );
	filters::Filters_map::const_iterator filter_it( filters.find( filter_name ) );
	if ( filter_it == filters.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Filter "+filter_name+" not found" );
	}
	filter_ =  filter_it->second;

	lower_bound_ = tag->getOption<Real>( "lower_bound");
	upper_bound_ = tag->getOption<Real>( "upper_bound");
	assert(lower_bound_ < upper_bound_);
}

filters::FilterOP
RangeFilterCreator::create_filter() const { return filters::FilterOP( new RangeFilter ); }

std::string
RangeFilterCreator::keyname() const { return "Range"; }


} // filters
} // protocols

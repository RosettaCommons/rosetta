// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/TimeFilter.cc
/// @author Sarel Fleishman (sarelf@uw.edu)

#include <protocols/filters/TimeFilter.hh>
#include <protocols/filters/TimeFilterCreator.hh>

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

//parsing
// AUTO-REMOVED #include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <time.h>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace filters {

static thread_local basic::Tracer TR( "protocols.filters.TimeFilter" );

///@brief default ctor
TimeFilter::TimeFilter() :
	Filter( "Time" ),
	tic_( 0 ),
	toc_( 0 )
{}

TimeFilter::~TimeFilter() {}

core::Real
TimeFilter::compute( core::pose::Pose const & ) const
{
	if( tic_ == 0 ){
		tic_ = time( NULL );
	}
	else{
		toc_ = time( NULL );
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

FilterOP
TimeFilterCreator::create_filter() const { return new TimeFilter; }

std::string
TimeFilterCreator::keyname() const { return "Time"; }

} // filters
} // protocols

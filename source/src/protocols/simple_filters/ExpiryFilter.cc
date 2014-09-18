// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ExpiryFilter.cc
/// @brief
/// @author Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/ExpiryFilter.hh>
#include <protocols/simple_filters/ExpiryFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <time.h>
namespace protocols{
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static thread_local basic::Tracer TR( "protocols.simple_filters.ExpiryFilter" );

protocols::filters::FilterOP
ExpiryFilterCreator::create_filter() const { return new ExpiryFilter; }

std::string
ExpiryFilterCreator::keyname() const { return "Expiry"; }

//default ctor
ExpiryFilter::ExpiryFilter() :
protocols::filters::Filter( "Expiry" ),
seconds_( 0 )
{
	start_time_ = time( NULL );
}

ExpiryFilter::~ExpiryFilter() {}

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
	return( time( NULL ) - start_time() );
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

}
}

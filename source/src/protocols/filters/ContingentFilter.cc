// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)

#include <protocols/filters/ContingentFilter.hh>
#include <protocols/filters/ContingentFilterCreator.hh>

#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <utility/tag/Tag.hh> // REQUIRED FOR WINDOWS
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace filters {

static thread_local basic::Tracer TR( "protocols.filters.ContingentFilter" );

///@brief default ctor
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

FilterOP
ContingentFilterCreator::create_filter() const { return new ContingentFilter; }

std::string
ContingentFilterCreator::keyname() const { return "ContingentFilter"; }



} // filters
} // protocols

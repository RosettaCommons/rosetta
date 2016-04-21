// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

--class--::--class--():
	protocols::filters::Filter( "--class--" )
{

}

--class--::~--class--()
{}

void
--class--::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}

protocols::filters::FilterOP
--class--::clone() const
{
	return protocols::filters::FilterOP( new --class--( *this ) );
}


protocols::filters::FilterOP
--class--::fresh_instance() const
{
	return protocols::filters::FilterOP( new --class-- );
}

std::string
--class--::get_name() const
{
	return --class--Creator::filter_name();
}

bool
--class--::apply( core::pose::Pose const & ) const
{
	return true;
}

core::Real
--class--::report_sm( core::pose::Pose const & ) const
{
	return -99999.9;
}

void
--class--::report( std::ostream &, core::pose::Pose const & ) const
{

}

/////////////// Creator ///////////////

protocols::filters::FilterOP
--class--Creator::create_filter() const
{
	return protocols::filters::FilterOP( new --class-- );
}

std::string
--class--Creator::keyname() const
{
	return --class--Creator::filter_name();
}

std::string
--class--Creator::filter_name()
{
	return "--class--";
}

--end_namespace--


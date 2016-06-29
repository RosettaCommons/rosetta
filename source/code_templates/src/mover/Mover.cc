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

// Unit headers
#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

--class--::--class--():
	protocols::moves::Mover( --class--::class_name() )
{

}

--class--::~--class--(){}

void
--class--::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}

protocols::moves::MoverOP
--class--::clone() const
{
	return protocols::moves::MoverOP( new --class--( *this ) );
}

protocols::moves::MoverOP
--class--::fresh_instance() const
{
	return protocols::moves::MoverOP( new --class-- );
}

std::string
--class--::get_name() const
{
	return --class--::class_name();
}

std::string
--class--::class_name()
{
	return "--class--";
}

void
--class--::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, --class-- const & mover )
{
	mover.show(os);
	return os;
}

void
--class--::apply( core::pose::Pose & )
{
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
--class--Creator::create_mover() const
{
	return protocols::moves::MoverOP( new --class-- );
}

std::string
--class--Creator::keyname() const
{
	return --class--::class_name();
}

--end_namespace--


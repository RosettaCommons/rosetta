// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/LoadUnboundRotMover.cc
/// @brief LoadUnboundRotMover methods implemented
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com

// Unit Headers
#include <protocols/simple_moves/LoadUnboundRotMover.hh>
#include <protocols/simple_moves/LoadUnboundRotMoverCreator.hh>

// Package Headers
#include <core/pack/dunbrack/RotamerConstraint.hh>

// Project Headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/Tracer.hh>

// C++ Headers

static thread_local basic::Tracer TR( "protocols.simple_moves.LoadUnboundRotMover" );

namespace protocols {
namespace simple_moves {

std::string
LoadUnboundRotMoverCreator::keyname() const
{
	return LoadUnboundRotMoverCreator::mover_name();
}

protocols::moves::MoverOP
LoadUnboundRotMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoadUnboundRotMover );
}

std::string
LoadUnboundRotMoverCreator::mover_name()
{ //this is the string name that was hardcoded into the Parser at SVN 46190
	return "load_unbound_rot";
}

LoadUnboundRotMover::LoadUnboundRotMover()
: protocols::moves::Mover("LoadUnboundRotMover")
{}

LoadUnboundRotMover::~LoadUnboundRotMover(){}

/// @details
void LoadUnboundRotMover::apply( core::pose::Pose & pose ){
	core::pack::dunbrack::load_unboundrot(pose);
	return;
}//apply

std::string
LoadUnboundRotMover::get_name() const {
	return LoadUnboundRotMoverCreator::mover_name();
}

protocols::moves::MoverOP LoadUnboundRotMover::fresh_instance() const { return protocols::moves::MoverOP( new LoadUnboundRotMover ); }
protocols::moves::MoverOP LoadUnboundRotMover::clone() const { return protocols::moves::MoverOP( new LoadUnboundRotMover( *this ) ); }

/// @brief parse XML (specifically in the context of the parser/scripting scheme); it's a no-op
void
LoadUnboundRotMover::parse_my_tag(
	utility::tag::TagCOP const,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{}

}//simple_moves
}//protocols


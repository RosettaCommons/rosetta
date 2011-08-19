// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SwitchResidueTypeSetMover.cc
/// @brief switch between residue type sets (e.g. centroid and all atom)

// Unit headers
#include <protocols/moves/SwitchResidueTypeSetMover.hh>
#include <protocols/moves/SwitchResidueTypeSetMoverCreator.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>


#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("protocols.moves.SwitchResidueTypeSetMover");
#include <utility/tag/Tag.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>


namespace protocols {
namespace moves {

std::string
SwitchResidueTypeSetMoverCreator::keyname() const
{
	return SwitchResidueTypeSetMoverCreator::mover_name();
}

protocols::moves::MoverOP
SwitchResidueTypeSetMoverCreator::create_mover() const {
	return new SwitchResidueTypeSetMover;
}

std::string
SwitchResidueTypeSetMoverCreator::mover_name()
{
	return "SwitchResidueTypeSetMover";
}

SwitchResidueTypeSetMover::SwitchResidueTypeSetMover()
	: Mover("SwitchResidueTypeSetMover")
{}

SwitchResidueTypeSetMover::SwitchResidueTypeSetMover( std::string const & type_set_tag_in )
	: Mover("SwitchResidueTypeSetMover"),
		type_set_tag_( type_set_tag_in )
{}

void
SwitchResidueTypeSetMover::apply( Pose & pose )
{
	core::util::switch_to_residue_type_set( pose, type_set_tag_ );
}

std::string
SwitchResidueTypeSetMover::get_name() const {
	return SwitchResidueTypeSetMoverCreator::mover_name();
}

MoverOP
SwitchResidueTypeSetMover::clone() const
{
	return new SwitchResidueTypeSetMover( *this );
}

MoverOP
SwitchResidueTypeSetMover::fresh_instance() const
{
	return new SwitchResidueTypeSetMover;
}

void
SwitchResidueTypeSetMover::parse_my_tag(
	TagPtr const tag,
	DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	if ( tag->hasOption("set") ) type_set_tag_ = tag->getOption<std::string>("set");
}

} // moves
} // protocols

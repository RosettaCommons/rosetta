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
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMoverCreator.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>


#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.simple_moves.SwitchResidueTypeSetMover" );
#include <utility/tag/Tag.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/Jump.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

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
	: moves::Mover("SwitchResidueTypeSetMover")
{}

SwitchResidueTypeSetMover::SwitchResidueTypeSetMover( std::string const & type_set_tag_in )
	: moves::Mover("SwitchResidueTypeSetMover"),
		type_set_tag_( type_set_tag_in )
{}

std::string
SwitchResidueTypeSetMover::get_residue_type_set() const {
	return type_set_tag_;
}

void
SwitchResidueTypeSetMover::apply( Pose & pose )
{
	core::util::switch_to_residue_type_set( pose, type_set_tag_ );
}

std::string
SwitchResidueTypeSetMover::get_name() const {
	return SwitchResidueTypeSetMoverCreator::mover_name();
}

void
SwitchResidueTypeSetMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Residue type set: " << get_residue_type_set() << std::endl;
}

moves::MoverOP
SwitchResidueTypeSetMover::clone() const
{
	return new SwitchResidueTypeSetMover( *this );
}

moves::MoverOP
SwitchResidueTypeSetMover::fresh_instance() const
{
	return new SwitchResidueTypeSetMover;
}

void
SwitchResidueTypeSetMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	if ( tag->hasOption("set") ) type_set_tag_ = tag->getOption<std::string>("set");
}

void
SwitchResidueTypeSetMover::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & /*score_fxns*/,
				utility::lua::LuaObject const & /*tasks*/,
				protocols::moves::MoverCacheSP /*cache*/ ) {
  if( def["set"] ) type_set_tag_ = def["set"].to<std::string>();
}

std::ostream &operator<< (std::ostream &os, SwitchResidueTypeSetMover const &mover)
{
	mover.show(os);
	return os;
}

} // simple_moves
} // protocols

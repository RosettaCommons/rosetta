// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MinMover.cc
/// @brief
/// @author Ingemar Andre

// Unit headers
#include <protocols/symmetric_docking/SymFoldandDockMoveRbJumpMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

// Package headers
#include <core/pose/symmetry/util.hh>

#include <protocols/symmetric_docking/SymFoldandDockCreators.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace symmetric_docking {

static basic::Tracer TR( "protocols.moves.symmetry.SymFoldandDockMoveRbJumpMover" );


SymFoldandDockMoveRbJumpMover::SymFoldandDockMoveRbJumpMover()
: Mover("SymFoldandDockMoveRbJumpMover")
{}

void
SymFoldandDockMoveRbJumpMover::apply( core::pose::Pose & pose )
{
	using namespace core::conformation::symmetry;
	protocols::simple_moves::symmetry::SetupForSymmetryMover setup;
	setup.apply( pose );
	core::pose::symmetry::find_new_symmetric_jump_residues( pose );
}

// XRW TEMP std::string
// XRW TEMP SymFoldandDockMoveRbJumpMover::get_name() const {
// XRW TEMP  return "SymFoldandDockMoveRbJumpMover";
// XRW TEMP }

void
SymFoldandDockMoveRbJumpMover::parse_my_tag(
	utility::tag::TagCOP const,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
}


// XRW TEMP std::string
// XRW TEMP SymFoldandDockMoveRbJumpMoverCreator::keyname() const {
// XRW TEMP  return SymFoldandDockMoveRbJumpMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SymFoldandDockMoveRbJumpMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SymFoldandDockMoveRbJumpMover() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SymFoldandDockMoveRbJumpMover::mover_name() {
// XRW TEMP  return "SymFoldandDockMoveRbJumpMover";
// XRW TEMP }

std::string SymFoldandDockMoveRbJumpMover::get_name() const {
	return mover_name();
}

std::string SymFoldandDockMoveRbJumpMover::mover_name() {
	return "SymFoldandDockMoveRbJumpMover";
}

void SymFoldandDockMoveRbJumpMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This seems to help set up Symmetric Fold and Dock. It symmetrizes the asymmetric input pose and then finds the jump residues between rigid bodies.", attlist );
}

std::string SymFoldandDockMoveRbJumpMoverCreator::keyname() const {
	return SymFoldandDockMoveRbJumpMover::mover_name();
}

protocols::moves::MoverOP
SymFoldandDockMoveRbJumpMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SymFoldandDockMoveRbJumpMover );
}

void SymFoldandDockMoveRbJumpMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymFoldandDockMoveRbJumpMover::provide_xml_schema( xsd );
}


} // symmetric_docking
} // protocols

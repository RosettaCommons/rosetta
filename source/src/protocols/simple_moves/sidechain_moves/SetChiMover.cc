// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/sidechain_moves/SetChiMover.cc
/// @brief  A mover to change one chi angle
/// @author Noah Ollikanen

// Unit headers
#include <protocols/simple_moves/sidechain_moves/SetChiMover.hh>
#include <protocols/simple_moves/sidechain_moves/SetChiMoverCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
//parsing
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

using namespace core;
using namespace core::chemical;
using namespace std;

using core::pose::Pose;
using core::conformation::Residue;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.sidechain_moves.SetChiMover" );

// XRW TEMP std::string
// XRW TEMP SetChiMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SetChiMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SetChiMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SetChiMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SetChiMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "SetChiMover";
// XRW TEMP }

SetChiMover::~SetChiMover() {}

/// @brief default ctor
SetChiMover::SetChiMover() :
	parent(),
	angle_( 0 ),
	resnum_( 0 ),
	chinum_( 0 )
{}

void SetChiMover::apply( Pose & pose ) {
	runtime_assert( resnum() > 0 );
	runtime_assert( resnum() <= pose.size() );

	if ( chinum() <= pose.residue(resnum()).nchi() ) {
		pose.set_chi(chinum(), resnum(), angle());
		TR<<"Set chi"<<chinum()<<" of residue "<<resnum()<<" to "<<angle()<<std::endl;
	}

	pose.update_residue_neighbors();
}

// XRW TEMP std::string
// XRW TEMP SetChiMover::get_name() const {
// XRW TEMP  return SetChiMover::mover_name();
// XRW TEMP }

void SetChiMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & pose)
{
	angle( tag->getOption< core::Real >( "angle" ) );
	resnum( core::pose::parse_resnum( tag->getOption< std::string >( "resnum" ), pose ) );
	chinum( tag->getOption< core::Size >( "chinum" ) );

}

std::string SetChiMover::get_name() const {
	return mover_name();
}

std::string SetChiMover::mover_name() {
	return "SetChiMover";
}

void SetChiMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// AMW no proper "defaults" but no requirement.
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute("angle", xsct_real, "Angle to which to set the chi." )
		+ XMLSchemaAttribute("resnum", xsct_refpose_enabled_residue_number, "Residue with the chi in question." )
		+ XMLSchemaAttribute("chinum", xsct_non_negative_integer, "Which chi is to be set." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "A mover to change one chi angle.", attlist );
}

std::string SetChiMoverCreator::keyname() const {
	return SetChiMover::mover_name();
}

protocols::moves::MoverOP
SetChiMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetChiMover );
}

void SetChiMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetChiMover::provide_xml_schema( xsd );
}



} // sidechain_moves
} // simple_moves
} // protocols

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Remove constraints from the current pose conformation.
/// @author Yifan Song

#include <protocols/simple_moves/ClearConstraintsMover.hh>
#include <protocols/simple_moves/ClearConstraintsMoverCreator.hh>

#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

ClearConstraintsMover::ClearConstraintsMover()= default;
ClearConstraintsMover::~ClearConstraintsMover()= default;

void ClearConstraintsMover::apply( core::pose::Pose & pose )
{
	pose.remove_constraints();
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
ClearConstraintsMover::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
}

moves::MoverOP ClearConstraintsMover::clone() const { return moves::MoverOP( new ClearConstraintsMover( *this ) ); }
moves::MoverOP ClearConstraintsMover::fresh_instance() const { return moves::MoverOP( new ClearConstraintsMover ); }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ClearConstraintsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ClearConstraintsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ClearConstraintsMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ClearConstraintsMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ClearConstraintsMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "ClearConstraintsMover";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ClearConstraintsMover::get_name() const {
// XRW TEMP  return "ClearConstraintsMover";
// XRW TEMP }

std::string ClearConstraintsMover::get_name() const {
	return mover_name();
}

std::string ClearConstraintsMover::mover_name() {
	return "ClearConstraintsMover";
}

void ClearConstraintsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Removes constraints from the Pose", attlist );
}

std::string ClearConstraintsMoverCreator::keyname() const {
	return ClearConstraintsMover::mover_name();
}

protocols::moves::MoverOP
ClearConstraintsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ClearConstraintsMover );
}

void ClearConstraintsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ClearConstraintsMover::provide_xml_schema( xsd );
}


} // moves
} // protocols

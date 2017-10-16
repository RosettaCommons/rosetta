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
/// @author Christoffer Norn

#include <protocols/simple_moves/SaveSequenceToCommentsMover.hh>
#include <protocols/simple_moves/SaveSequenceToCommentsMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <utility/tag/Tag.hh>


// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

SaveSequenceToCommentsMover::SaveSequenceToCommentsMover() :
	save_seq_name_("saved_seq")
{
}

SaveSequenceToCommentsMover::~SaveSequenceToCommentsMover()= default;

void SaveSequenceToCommentsMover::apply( core::pose::Pose & pose )
{
	std::string aa_seq = pose.sequence();
	core::pose::add_comment(pose, save_seq_name(), aa_seq);
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
SaveSequenceToCommentsMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
	save_seq_name( tag->getOption< std::string >( "save_seq_name", "" ) );
}

moves::MoverOP SaveSequenceToCommentsMover::clone() const { return moves::MoverOP( new SaveSequenceToCommentsMover( *this ) ); }
moves::MoverOP SaveSequenceToCommentsMover::fresh_instance() const { return moves::MoverOP( new SaveSequenceToCommentsMover ); }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SaveSequenceToCommentsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SaveSequenceToCommentsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SaveSequenceToCommentsMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SaveSequenceToCommentsMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SaveSequenceToCommentsMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "SaveSequenceToCommentsMover";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SaveSequenceToCommentsMover::get_name() const {
// XRW TEMP  return "SaveSequenceToCommentsMover";
// XRW TEMP }

std::string SaveSequenceToCommentsMover::get_name() const {
	return mover_name();
}

std::string SaveSequenceToCommentsMover::mover_name() {
	return "SaveSequenceToCommentsMover";
}

void SaveSequenceToCommentsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "save_seq_name", xs_string, "Name under which to save the current save to comments" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Saves the current sequence to comments", attlist );
}

std::string SaveSequenceToCommentsMoverCreator::keyname() const {
	return SaveSequenceToCommentsMover::mover_name();
}

protocols::moves::MoverOP
SaveSequenceToCommentsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SaveSequenceToCommentsMover );
}

void SaveSequenceToCommentsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SaveSequenceToCommentsMover::provide_xml_schema( xsd );
}


} // moves
} // protocols

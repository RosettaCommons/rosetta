// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DeleteChainMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/DeleteChainMover.hh>
#include <protocols/simple_moves/DeleteChainMoverCreator.hh>

// Core headers
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.DeleteChainMover" );

// XRW TEMP std::string
// XRW TEMP DeleteChainMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DeleteChainMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DeleteChainMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new DeleteChainMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DeleteChainMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "DeleteChain";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DeleteChainMover::get_name() const {
// XRW TEMP  return DeleteChainMover::mover_name();
// XRW TEMP }

moves::MoverOP
DeleteChainMover::clone() const
{
	return moves::MoverOP( new DeleteChainMover( *this ) );
}

moves::MoverOP
DeleteChainMover::fresh_instance() const
{
	return moves::MoverOP( new DeleteChainMover );
}

DeleteChainMover::DeleteChainMover():
	moves::Mover("DeleteChain"),
	chain_num_( 1 )
{}

///////////////END BOILER PLATE CODE////////////////

void
DeleteChainMover::chain_num(
	core::Size chain_num
) {
	chain_num_ = chain_num;
}

core::Size
DeleteChainMover::chain_num() {
	return chain_num_;
}

void
DeleteChainMover::apply( Pose & pose )
{
	utility::vector1<core::pose::PoseOP> chain_poses = pose.split_by_chain();
	TR << "Removing chain " << chain_num() << " from pose with " << chain_poses.size() << " chains" << std::endl;

	bool chain_added = false;
	for ( core::Size i=1; i<=chain_poses.size(); ++i ) {
		if ( i != chain_num_ ) {
			if ( !chain_added ) {
				pose = *chain_poses[i];
				chain_added = true;
			} else {
				pose.append_pose_by_jump(*chain_poses[i], 1);
			}
		}
	}
}


void
DeleteChainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	if ( tag->hasOption("chain") ) {
		chain_num_ = tag->getOption< core::Size >( "chain");
	}
}

std::string DeleteChainMover::get_name() const {
	return mover_name();
}

std::string DeleteChainMover::mover_name() {
	return "DeleteChain";
}

void DeleteChainMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute( "chain", xsct_non_negative_integer, "delete this chain (number) (no option for PDB chains)");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Delete this chain.", attlist );
}

std::string DeleteChainMoverCreator::keyname() const {
	return DeleteChainMover::mover_name();
}

protocols::moves::MoverOP
DeleteChainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DeleteChainMover );
}

void DeleteChainMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DeleteChainMover::provide_xml_schema( xsd );
}


} // simple_moves
} // protocols

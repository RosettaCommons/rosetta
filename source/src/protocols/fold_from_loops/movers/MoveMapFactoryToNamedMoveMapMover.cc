// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/MoveMapFactoryToNamedMoveMapMover.cc
/// @brief Adds to the DataMap a MoveMap obtained from applying the MoveMapFactory
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/movers/MoveMapFactoryToNamedMoveMapMover.hh>
#include <protocols/fold_from_loops/movers/MoveMapFactoryToNamedMoveMapMoverCreator.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/util.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

static basic::Tracer TR( "protocols.fold_from_loops.MoveMapFactoryToNamedMoveMapMover", basic::t_trace );

namespace protocols {
namespace fold_from_loops {
namespace movers {

MoveMapFactoryToNamedMoveMapMover::MoveMapFactoryToNamedMoveMapMover():
	protocols::moves::Mover( mover_name() ),
	movemap_( new core::kinematics::MoveMap )
{}

MoveMapFactoryToNamedMoveMapMover::~MoveMapFactoryToNamedMoveMapMover()= default;

void
MoveMapFactoryToNamedMoveMapMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	movemapfactory_ = core::select::movemap::parse_movemap_factory( tag, data );
	std::string const movemapname( tag->getOption< std::string >( "movemap" ) );
	if ( !data.has("movemaps", movemapname) ) {
		data.add( "movemaps", movemapname, movemap_ );
	}
}

protocols::moves::MoverOP
MoveMapFactoryToNamedMoveMapMover::clone() const
{
	return protocols::moves::MoverOP( new MoveMapFactoryToNamedMoveMapMover( *this ) );
}

protocols::moves::MoverOP
MoveMapFactoryToNamedMoveMapMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new MoveMapFactoryToNamedMoveMapMover );
}

void
MoveMapFactoryToNamedMoveMapMover::apply( core::pose::Pose & pose )
{
	core::kinematics::MoveMapOP tmpmap = movemapfactory_->create_movemap_from_pose( pose );
	*movemap_ = *tmpmap;

}

std::string MoveMapFactoryToNamedMoveMapMover::get_name() const {
	return mover_name();
}

std::string MoveMapFactoryToNamedMoveMapMover::mover_name() {
	return "MoveMapFactoryToNamedMoveMapMover";
}

void MoveMapFactoryToNamedMoveMapMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "movemap", xs_string, "Name of the expected created MoveMap");
	core::select::movemap::attributes_for_parse_movemap_factory_when_required( attlist, "movemap_factory", "Name of the MoveMapFactory to be Processed");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Adds to the DataMap a MoveMap obtained from applying the MoveMapFactory.", attlist );
}

std::string MoveMapFactoryToNamedMoveMapMoverCreator::keyname() const {
	return MoveMapFactoryToNamedMoveMapMover::mover_name();
}

protocols::moves::MoverOP
MoveMapFactoryToNamedMoveMapMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MoveMapFactoryToNamedMoveMapMover );
}

void MoveMapFactoryToNamedMoveMapMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MoveMapFactoryToNamedMoveMapMover::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops

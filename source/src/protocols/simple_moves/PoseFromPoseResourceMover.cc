// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/PoseFromPoseResourceMover.cc
/// @brief Mover that shuttles a Pose from a PoseResource into the DataMap when its parse_my_tag method is invoked
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/simple_moves/PoseFromPoseResourceMover.hh>
#include <protocols/simple_moves/PoseFromPoseResourceMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/PoseResource.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.simple_moves.PoseFromPoseResourceMover" );

namespace protocols {
namespace simple_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PoseFromPoseResourceMover::PoseFromPoseResourceMover():
	protocols::moves::Mover( PoseFromPoseResourceMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PoseFromPoseResourceMover::PoseFromPoseResourceMover( PoseFromPoseResourceMover const & src ):
	protocols::moves::Mover( src )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PoseFromPoseResourceMover::~PoseFromPoseResourceMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

void
PoseFromPoseResourceMover::apply( core::pose::Pose& )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
PoseFromPoseResourceMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
PoseFromPoseResourceMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	std::string pose_resource_name = tag->getOption< std::string >( "pose_resource_name" );
	core::import_pose::PoseResourceCOP resource;
	try {
		resource = data_map.get_resource< core::import_pose::PoseResource >( pose_resource_name );
	} catch ( utility::excn::Exception & e ) {
		std::ostringstream oss;
		oss << "Failed to retrieve PoseResource from the DataMap in PoseFromPoseResourceMover::parse_my_tag\n";
		oss << "Error message below:\n" << e.msg();
		throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	TR << "Adding a copy of the PoseResource named '" << pose_resource_name << "' to the data map" << std::endl;
	data_map.add( "poses", pose_resource_name, resource->pose_deep_copy() );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PoseFromPoseResourceMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PoseFromPoseResourceMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PoseFromPoseResourceMover::clone() const
{
	return protocols::moves::MoverOP( new PoseFromPoseResourceMover( *this ) );
}

std::string PoseFromPoseResourceMover::get_name() const {
	return mover_name();
}

std::string PoseFromPoseResourceMover::mover_name() {
	return "PoseFromPoseResourceMover";
}

void PoseFromPoseResourceMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute( "pose_resource_name", xs_string, "The name of the PoseResource"
		"object that must be retrieved from the DataMap's resources map." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "During its creation within the rosetta_scripts"
		" application (specifically, during the invocation of its parse_my_tag method), this mover will extract a PoseResource"
		" from the DataMap's map of (bitwise-const) Resources and from that PoseResource, extract a deep copy of the Pose"
		" that it holds. It will then place that Pose into the DataMap so that a single copy of that Pose can be shared by"
		" many Mover/Filters. Why the dance? Why not just hold on to the PoseResource itself? The complication with Pose is that"
		" it is not a bitwise-const class -- that is, it holds mutable data. If a single Pose were shared by two threads, then"
		" race conditions would likely ensue. To fit into the definition of a bitwise-const object,"
		" the PoseResource only allows you to access a deep copy of the original Pose that it holds. If the Pose is expensive to"
		" construct, then loading it in from disk a single time might be advantageous. But then, you perhaps do not want to have"
		" a separate copy of that Pose for each object that is going to want to access that Pose. The purpose of this mover is"
		" to shuttle a Pose out of the safe-to-be-shared-between-threads Resource section of the DataMap into the only-accessed-"
		"by-a-single-thread section of the DataMap.", attlist );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
PoseFromPoseResourceMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new PoseFromPoseResourceMover );
}

std::string
PoseFromPoseResourceMoverCreator::keyname() const
{
	return PoseFromPoseResourceMover::mover_name();
}

void PoseFromPoseResourceMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PoseFromPoseResourceMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, PoseFromPoseResourceMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //simple_moves

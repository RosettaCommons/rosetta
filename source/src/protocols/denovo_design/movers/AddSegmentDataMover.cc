// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/AddSegmentDataMover.cc
/// @brief Adds a segment to the structuredata
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/movers/AddSegmentDataMover.hh>
#include <protocols/denovo_design/movers/AddSegmentDataMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.AddSegmentDataMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

AddSegmentDataMover::AddSegmentDataMover():
	protocols::moves::Mover( AddSegmentDataMover::mover_name() ),
	segment_name_( "" ),
	secstruct_( "" ),
	abego_( "" )
{

}

AddSegmentDataMover::~AddSegmentDataMover(){}

void
AddSegmentDataMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	segment_name_ = tag->getOption< std::string >( "segment_name", segment_name_ );
	abego_ = tag->getOption< std::string >( "abego", abego_ );
	secstruct_ = tag->getOption< std::string >( "secstruct", secstruct_ );

	// input checking
	if ( segment_name_.empty() ) {
		std::stringstream msg;
		msg << "Name of new segment must be specified to AddSegment" << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	if ( secstruct_.empty() || abego_.empty() ) {
		std::stringstream msg;
		msg << "SS and ABEGO must be specified to AddSegment" << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	if ( secstruct_.size() != abego_.size() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "AddSegment: ss and abego must be the same length!" );
	}
}

protocols::moves::MoverOP
AddSegmentDataMover::clone() const
{
	return protocols::moves::MoverOP( new AddSegmentDataMover( *this ) );
}

protocols::moves::MoverOP
AddSegmentDataMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AddSegmentDataMover );
}

// XRW TEMP std::string
// XRW TEMP AddSegmentDataMover::get_name() const
// XRW TEMP {
// XRW TEMP  return AddSegmentDataMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddSegmentDataMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "AddSegmentDataMover";
// XRW TEMP }

void
AddSegmentDataMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, AddSegmentDataMover const & mover )
{
	mover.show(os);
	return os;
}

void
AddSegmentDataMover::apply( core::pose::Pose & pose )
{
	components::StructureData newsd = components::StructureDataFactory::get_instance()->get_from_pose( pose );
	create_segment( newsd );
	components::StructureDataFactory::get_instance()->save_into_pose( pose, newsd );
}

void
AddSegmentDataMover::create_segment( components::StructureData & perm ) const
{
	if ( segment_name_.empty() ) {
		std::stringstream msg;
		msg << "Name of new segment must be specified to AddSegment" << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	if ( secstruct_.empty() || abego_.empty() ) {
		std::stringstream msg;
		msg << "SS and ABEGO must be specified to AddSegment" << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}

	if ( secstruct_.size() != abego_.size() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "AddSegment: ss and abego must be the same length!" );
	}

	bool const nterm_inc = secstruct_.size() < 2;
	bool const cterm_inc = secstruct_.size() < 3;
	components::Segment const seg( segment_name_, secstruct_, abego_, nterm_inc, cterm_inc );
	perm.add_segment( seg );
}

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddSegmentDataMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new AddSegmentDataMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddSegmentDataMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AddSegmentDataMover::mover_name();
// XRW TEMP }

std::string AddSegmentDataMover::get_name() const {
	return mover_name();
}

std::string AddSegmentDataMover::mover_name() {
	return "AddSegmentDataMover";
}

void AddSegmentDataMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute( "segment_name", xs_string, "Name of new segment")
		+ XMLSchemaAttribute::required_attribute( "abego", xs_string, "ABEGO of new segment")
		+ XMLSchemaAttribute::required_attribute( "secstruct", xsct_dssp_string, "Secondary structure of new segment. Must be same length as ABEGO string."); //there is a dependency between abego and ss (must be same length)
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Adds a segment to the StructureData", attlist );

}

std::string AddSegmentDataMoverCreator::keyname() const {
	return AddSegmentDataMover::mover_name();
}

protocols::moves::MoverOP
AddSegmentDataMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddSegmentDataMover );
}

void AddSegmentDataMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddSegmentDataMover::provide_xml_schema( xsd );
}


} //protocols
} //denovo_design
} //movers


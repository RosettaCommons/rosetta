// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.AddSegmentDataMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

AddSegmentDataMover::AddSegmentDataMover():
	protocols::moves::Mover( AddSegmentDataMover::class_name() ),
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

std::string
AddSegmentDataMover::get_name() const
{
	return AddSegmentDataMover::class_name();
}

std::string
AddSegmentDataMover::class_name()
{
	return "AddSegmentData";
}

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
	components::Segment const seg( secstruct_, abego_, nterm_inc, cterm_inc );
	perm.add_segment( segment_name_, seg );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
AddSegmentDataMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AddSegmentDataMover );
}

std::string
AddSegmentDataMoverCreator::keyname() const
{
	return AddSegmentDataMover::class_name();
}

} //protocols
} //denovo_design
} //movers


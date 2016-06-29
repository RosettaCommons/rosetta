// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/movers/SetResidueAliasMover.cc
/// @brief Sets a residue alias in the StructureData
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/movers/SetResidueAliasMover.hh>
#include <protocols/denovo_design/movers/SetResidueAliasMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.SetResidueAliasMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

SetResidueAliasMover::SetResidueAliasMover():
	protocols::moves::Mover( SetResidueAliasMover::class_name() ),
	alias_name_( "" ),
	segment_name_( "" ),
	resid_( 0 )
{

}

SetResidueAliasMover::~SetResidueAliasMover(){}

void
SetResidueAliasMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	alias_name_ = tag->getOption< std::string >( "alias_name", alias_name_ );
	segment_name_ = tag->getOption< std::string >( "segment_name", segment_name_ );
	resid_ = tag->getOption< core::Size >( "residue", resid_ );

	if ( resid_ == 0 ) {
		std::stringstream msg;
		msg << class_name() << "::apply(): Residue number must be specified to add an alias via the \"residue\" option. If a residue number is "
			<< "specified along with a segment_name, the residue number will be local within the segment." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}

	if ( alias_name_.empty() ) {
		std::stringstream msg;
		msg << class_name() << "::apply(): Alias name is not specified -- you must specify the \"alias_name\" option in order "
			<< "to set an alias." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
}

protocols::moves::MoverOP
SetResidueAliasMover::clone() const
{
	return protocols::moves::MoverOP( new SetResidueAliasMover( *this ) );
}

protocols::moves::MoverOP
SetResidueAliasMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new SetResidueAliasMover );
}

std::string
SetResidueAliasMover::get_name() const
{
	return SetResidueAliasMover::class_name();
}

std::string
SetResidueAliasMover::class_name()
{
	return "SetResidueAlias";
}

void
SetResidueAliasMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, SetResidueAliasMover const & mover )
{
	mover.show(os);
	return os;
}

void
SetResidueAliasMover::apply( core::pose::Pose & pose )
{
	if ( resid_ == 0 ) {
		std::stringstream msg;
		msg << class_name() << "::apply(): Residue number must be specified to add an alias. If a residue number is "
			<< "specified along with a segment_name, the residue number will be local within the segment." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( alias_name_.empty() ) {
		std::stringstream msg;
		msg << class_name() << "::apply(): Alias name is not specified -- you must specify the alias_name option in order "
			<< "to set an alias." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	components::StructureData sd = components::StructureDataFactory::get_instance()->get_from_pose( pose );
	if ( segment_name_.empty() ) {
		sd.set_alias( alias_name_, resid_ );
	} else {
		sd.set_alias( alias_name_, segment_name_, resid_ );
	}

	components::StructureDataFactory::get_instance()->save_into_pose( pose, sd );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
SetResidueAliasMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new SetResidueAliasMover );
}

std::string
SetResidueAliasMoverCreator::keyname() const
{
	return SetResidueAliasMover::class_name();
}

} //protocols
} //denovo_design
} //movers


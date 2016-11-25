// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.SetResidueAliasMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

SetResidueAliasMover::SetResidueAliasMover():
	protocols::moves::Mover( SetResidueAliasMover::mover_name() ),
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
		msg << mover_name() << "::apply(): Residue number must be specified to add an alias via the \"residue\" option. If a residue number is "
			<< "specified along with a segment_name, the residue number will be local within the segment." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}

	if ( alias_name_.empty() ) {
		std::stringstream msg;
		msg << mover_name() << "::apply(): Alias name is not specified -- you must specify the \"alias_name\" option in order "
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

// XRW TEMP std::string
// XRW TEMP SetResidueAliasMover::get_name() const
// XRW TEMP {
// XRW TEMP  return SetResidueAliasMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SetResidueAliasMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "SetResidueAlias";
// XRW TEMP }

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
		msg << mover_name() << "::apply(): Residue number must be specified to add an alias. If a residue number is "
			<< "specified along with a segment_name, the residue number will be local within the segment." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( alias_name_.empty() ) {
		std::stringstream msg;
		msg << mover_name() << "::apply(): Alias name is not specified -- you must specify the alias_name option in order "
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

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SetResidueAliasMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new SetResidueAliasMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SetResidueAliasMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SetResidueAliasMover::mover_name();
// XRW TEMP }

std::string SetResidueAliasMover::get_name() const {
	return mover_name();
}

std::string SetResidueAliasMover::mover_name() {
	return "SetResidueAlias";
}

void SetResidueAliasMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist; // TO DO: add attributes to this list
	attlist + XMLSchemaAttribute::required_attribute( "alias_name", xs_string, "Alias name to set for the specified residue")
		+ XMLSchemaAttribute( "segment_name", xs_string, "Optional argument to specify which segment this residue is in. If specified, residue numbering is local within the segment.")
		+ XMLSchemaAttribute::required_attribute( "residue", xsct_non_negative_integer, "Number of residue for which to specify an alias");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Used to specify an alias name in the StructureDefinition for a specified residue", attlist );
}

std::string SetResidueAliasMoverCreator::keyname() const {
	return SetResidueAliasMover::mover_name();
}

protocols::moves::MoverOP
SetResidueAliasMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetResidueAliasMover );
}

void SetResidueAliasMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetResidueAliasMover::provide_xml_schema( xsd );
}


} //protocols
} //denovo_design
} //movers


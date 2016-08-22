// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/CompoundArchitect.cc
/// @brief Architect that creates a StructureData using multiple architects
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/CompoundArchitectCreator.hh>

// Protocol headers
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/connection/util.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.CompoundArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

CompoundArchitect::CompoundArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	architects_(),
	connections_()
{
}

CompoundArchitect::~CompoundArchitect()
{}

CompoundArchitect::DeNovoArchitectOP
CompoundArchitect::clone() const
{
	return DeNovoArchitectOP( new CompoundArchitect( *this ) );
}

std::string
CompoundArchitect::type() const
{
	return CompoundArchitect::class_name();
}

void
CompoundArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	for ( utility::tag::Tag::tags_t::const_iterator t=tag->getTags().begin(); t!=tag->getTags().end(); ++t ) {
		if ( (*t)->getName() == "Architects" ) parse_architect_tags( *t, data );
		else if ( (*t)->getName() == "Connections" ) parse_connection_tags( *t, data );
		else if ( (*t)->getName() == "Pairing" ) parse_pairing_tags( **t );
		else {
			std::stringstream msg;
			msg << type() << "::parse_tag(): The provided subtag name (" << (*t)->getName()
				<< ") is not valid. Valid subtag names for " << type() << " are: "
				<< "Architects, Connections, Pairing." << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
	}
}

void
CompoundArchitect::parse_architect_tags( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	for ( utility::tag::Tag::tags_t::const_iterator t=tag->getTags().begin(); t!=tag->getTags().end(); ++t ) {
		parse_architect_tag( *t, data );
	}
}

void
CompoundArchitect::parse_architect_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	DeNovoArchitectOP new_architect = DeNovoArchitectFactory::get_instance()->create_from_tag( tag, data );
	architects_.push_back( new_architect );
}

void
CompoundArchitect::parse_connection_tags( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	for ( utility::tag::Tag::tags_t::const_iterator t=tag->getTags().begin(); t!=tag->getTags().end(); ++t ) {
		parse_connection_tag( *t, data );
	}
}

void
CompoundArchitect::parse_connection_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	connection::ConnectionArchitectOP architect;
	if ( tag->getName() == "Connection" ) {
		architect = connection::ConnectionArchitectOP(
			new connection::ConnectionArchitect( tag->getOption< std::string >( "name" ) ) );
	} else {
		std::stringstream msg;
		msg << "Bad connection architect type in tag: " << *tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	architect->parse_my_tag( tag, data );
	if ( ! architect->id().empty() ) connection::store_connection_architect( architect, data );
	connections_.push_back( architect );
}

void
CompoundArchitect::parse_pairing_tags( utility::tag::Tag const & tag )
{
	for ( utility::tag::Tag::tags_t::const_iterator t=tag.getTags().begin(); t!=tag.getTags().end(); ++t ) {
		parse_pairing_tag( **t );
	}
}

void
CompoundArchitect::parse_pairing_tag( utility::tag::Tag const & tag )
{
	components::SegmentPairingOP pair = components::SegmentPairing::create( tag.getName() );
	pair->parse_my_tag( tag );
	pairings_.push_back( pair );
}

CompoundArchitect::StructureDataOP
CompoundArchitect::design( core::pose::Pose const & pose, core::Real & random ) const
{
	StructureDataOP sd( new StructureData( id() ) );
	for ( DeNovoArchitectCOPs::const_iterator a=architects_.begin(); a!=architects_.end(); ++a ) {
		StructureDataOP sub_sd = (*a)->design( pose, random );
		if ( !sub_sd ) {
			std::stringstream msg;
			msg << class_name() << "::design(): Architect " << (*a)->id()
				<< " did not produce a StructureData object. Current SD ="
				<< *sd << std::endl;
			utility_exit_with_message( msg.str() );
		}
		if ( sub_sd->segments_begin() == sub_sd->segments_end() ) {
			std::stringstream msg;
			msg << class_name() << "::design(): Architect " << (*a)->id()
				<< " produced a StructureData object without segments. Current SD ="
				<< *sd << " Architect-produced SD = " << *sub_sd << std::endl;
			utility_exit_with_message( msg.str() );
		}
		sd->merge( *sub_sd );
	}

	for ( connection::ConnectionArchitectCOPs::const_iterator c=connections_.begin(); c!=connections_.end(); ++c ) {
		(*c)->apply( *sd, random );
	}

	for ( components::SegmentPairingCOPs::const_iterator p=pairings_.begin(); p!=pairings_.end(); ++p ) {
		sd->add_pairing( **p );
	}
	return sd;
}

void
CompoundArchitect::add_architect( DeNovoArchitect const & architect )
{
	architects_.push_back( architect.clone() );
}

void
CompoundArchitect::add_connection( connection::ConnectionArchitect const & connection )
{
	connections_.push_back( connection.clone() );
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
std::string
CompoundArchitectCreator::keyname() const
{
	return CompoundArchitect::class_name();
}

DeNovoArchitectOP
CompoundArchitectCreator::create_architect( std::string const & architect_id ) const
{
	return DeNovoArchitectOP( new CompoundArchitect( architect_id ) );
}

} //protocols
} //denovo_design
} //architects

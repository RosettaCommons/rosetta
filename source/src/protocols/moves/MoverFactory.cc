// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MoverFactory.cc
/// @brief
/// @author ashworth

#include <protocols/moves/MoverFactory.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/mover_schemas.hh>

// required for passing to Mover::parse_my_tag

#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, throw utility::excn::EXCN_RosettaScriptsOption
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>


namespace protocols {
namespace moves {


static basic::Tracer TR( "protocols.moves.MoverFactory" );

MoverFactory::MoverFactory()
{
	forbidden_names_.clear();
	//the following 5 names are used in the original
	//LoopMoverFactory
	forbidden_names_.insert("quick_ccd");
	forbidden_names_.insert("sdwindow");
	forbidden_names_.insert("quick_ccd_moves");
	forbidden_names_.insert("perturb_ccd");
	forbidden_names_.insert("perturb_kic");
	//original loop mover names over
}

MoverFactory::~MoverFactory()= default;

/// @brief add a Mover prototype, using its default type name as the map key
void
MoverFactory::factory_register( MoverCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	std::string const mover_type( creator->keyname() );
	if ( mover_type == "UNDEFINED NAME" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Can't map derived Mover with undefined type name.");
	}
	if ( forbidden_names_.find( mover_type ) != forbidden_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Name "+mover_type+" is not an allowed mover name, probably because it has historical meaning.");
	}
	if ( mover_creator_map_.find( mover_type ) != mover_creator_map_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "MoverFactory::factory_register already has a mover creator with name \"" + mover_type + "\".  Conflicting Mover names" );
	}
	mover_creator_map_[ mover_type ] = creator;
}

/// @brief Is there a mover with the given name that's known to Rosetta?
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
MoverFactory::mover_exists(
	std::string const &mover_name
) const {
	return( mover_creator_map_.find( mover_name ) != mover_creator_map_.end() );
}

/// @brief Get the XML schema for a given mover.
/// @details Throws an error if the mover is unknown to Rosetta.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
MoverFactory::provide_xml_schema(
	std::string const &mover_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	auto iter( mover_creator_map_.find( mover_name ) );
	if ( iter != mover_creator_map_.end() ) {
		if ( ! iter->second ) {
			utility_exit_with_message( "Error: MoverCreatorOP prototype for " + mover_name + " is NULL!" );
		}
		iter->second->provide_xml_schema( xsd );
	} else {
		TR << "Available movers: ";
		for ( auto const & filt_it : mover_creator_map_ ) {
			TR << filt_it.first << ", ";
		}
		TR << std::endl;
		utility_exit_with_message( mover_name + " is not known to the MoverFactory. Was it registered via a MoverRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
	}
}

/// @brief return new Mover by key lookup in mover_prototype_map_ (new Mover parses Tag if provided)
MoverOP
MoverFactory::newMover( std::string const & mover_type )
{
	MoverMap::const_iterator iter( mover_creator_map_.find( mover_type ) );
	if ( iter != mover_creator_map_.end() ) {
		if ( ! iter->second ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Error: MoverCreatorOP prototype for " + mover_type + " is NULL!" );
		}
		// use of cloning method would be faithful to pre-initialized prototypes
		//return iter->second->clone();
		// fresh_instance prevents propagation of pre-initialized prototypes, which may be safer(?)
		return iter->second->create_mover();
	} else {
		TR << "Available movers: ";
		for ( MoverMap::const_iterator mover_it = mover_creator_map_.begin(); mover_it != mover_creator_map_.end(); ++mover_it ) {
			TR << mover_it->first << ", ";
		}
		TR << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  mover_type + " is not known to the MoverFactory. Was it registered via a MoverRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return nullptr;
	}
}

/// @brief return new Mover by Tag parsing
MoverOP
MoverFactory::newMover(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose )
{
	MoverOP mover( newMover( tag->getName() ) );
	runtime_assert( mover != nullptr );
	mover->parse_my_tag( tag, data, filters, movers, pose );
	return mover;
}

MoverFactory::MoverMap const & MoverFactory::mover_creator_map() const { return mover_creator_map_; }

void MoverFactory::define_mover_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			mover_creator_map_,
			mover_xml_schema_group_name(),
			& complex_type_name_for_mover,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for Movers from MoverFactory; offending class"
			" must call protocols::moves::complex_type_name_for_mover when defining"
			" its XML Schema\n" + e.msg() );
	}

}

std::string MoverFactory::mover_xml_schema_group_name()
{
	return "mover";
}

} //namespace moves
} //namespace protocols

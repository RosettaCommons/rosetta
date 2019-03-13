// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/enzymatic_movers/GlycosyltransferaseMover.cc
/// @brief   Method definitions for GlycosyltransferaseMover.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Unit headers
#include <protocols/enzymatic_movers/GlycosyltransferaseMover.hh>
#include <protocols/enzymatic_movers/GlycosyltransferaseMoverCreator.hh>
#include <protocols/enzymatic_movers/EnzymaticMover.hh>

// Project headers
#include <core/pose/carbohydrates/util.hh>

#include <protocols/moves/mover_schemas.hh>
//#include <protocols/rosetta_scripts/util.hh>

// Utility headers
//#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Basic header
//#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>


// Construct tracers.
static basic::Tracer TR( "protocols.enzymatic_movers.GlycosyltransferaseMover" );


namespace protocols {
namespace enzymatic_movers {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
GlycosyltransferaseMover::GlycosyltransferaseMover(): EnzymaticMover( "glycosyltransferases" )
{
	init();
}

// Copy constructor
GlycosyltransferaseMover::GlycosyltransferaseMover( GlycosyltransferaseMover const & object_to_copy ) :
	EnzymaticMover( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

// Assignment operator
GlycosyltransferaseMover &
GlycosyltransferaseMover::operator=( GlycosyltransferaseMover const & object_to_copy )
{
	// Abort self-assignment.
	if ( this != &object_to_copy ) {
		EnzymaticMover::operator=( object_to_copy );
		copy_data( *this, object_to_copy );
	}
	return *this;
}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
GlycosyltransferaseMover::register_options()
{
	//using namespace basic::options;

	//option.add_relevant( OptionKeys::foo::bar );

	EnzymaticMover::register_options();
}


// Mover methods
std::string
GlycosyltransferaseMover::get_name() const {
	return mover_name();
}

moves::MoverOP
GlycosyltransferaseMover::clone() const
{
	return utility::pointer::make_shared< GlycosyltransferaseMover >( *this );
}

moves::MoverOP
GlycosyltransferaseMover::fresh_instance() const
{
	return utility::pointer::make_shared< GlycosyltransferaseMover >();
}


void
GlycosyltransferaseMover::parse_my_tag(
	TagCOP /*tag*/,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/ )
{}

void
GlycosyltransferaseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Enzymatic mover to glycosylate a pose", attlist );
}

// Protected methods //////////////////////////////////////////////////////////
void
GlycosyltransferaseMover::perform_reaction(
	core::pose::Pose & input_pose,
	core::uint const sepos,
	std::string const & cosubstrate )
{
	core::pose::carbohydrates::glycosylate_pose(
		input_pose,
		get_reactive_site_sequence_position( sepos ),
		get_reactive_site_atom_name( sepos ),
		cosubstrate );
}


// Private methods ////////////////////////////////////////////////////////////
// Set command-line options.  (Called by init())
void
GlycosyltransferaseMover::set_commandline_options()
{
	//using namespace basic::options;
}

// Initialize data members from arguments.
void
GlycosyltransferaseMover::init()
{
	type( "GlycosyltransferaseMover" );
	set_commandline_options();

	// Set defaults.
	set_species( "h_sapiens" );
	set_enzyme( "generic_N-glycosyltransferase" );
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
GlycosyltransferaseMover::copy_data(
	GlycosyltransferaseMover & /*object_to_copy_to*/,
	GlycosyltransferaseMover const & /*object_to_copy_from*/ )
{}


// Creator methods ////////////////////////////////////////////////////////////
std::string
GlycosyltransferaseMoverCreator::keyname() const {
	return GlycosyltransferaseMover::mover_name();
}

// Return an up-casted owning pointer (MoverOP) to the mover.
protocols::moves::MoverOP
GlycosyltransferaseMoverCreator::create_mover() const {
	return utility::pointer::make_shared< GlycosyltransferaseMover >();
}

void
GlycosyltransferaseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycosyltransferaseMover::provide_xml_schema( xsd );
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that GlycosyltransferaseMover can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, GlycosyltransferaseMover const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace enzymatic_movers
}  // namespace protocols

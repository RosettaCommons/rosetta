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

// Utility headers
#include <utility/tag/XMLSchemaGeneration.hh>


namespace protocols {
namespace enzymatic_movers {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
GlycosyltransferaseMover::GlycosyltransferaseMover(): EnzymaticMover( "glycosyltransferases" )
{
	type( "GlycosyltransferaseMover" );
}

// Copy constructor
GlycosyltransferaseMover::GlycosyltransferaseMover( GlycosyltransferaseMover const & object_to_copy ) :
	EnzymaticMover( object_to_copy )
{}

// Assignment operator
GlycosyltransferaseMover &
GlycosyltransferaseMover::operator=( GlycosyltransferaseMover const & object_to_copy )
{
	// Abort self-assignment.
	if ( this != &object_to_copy ) {
		EnzymaticMover::operator=( object_to_copy );
	}
	return *this;
}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
GlycosyltransferaseMover::register_options()
{
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
GlycosyltransferaseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	EnzymaticMover::xml_schema_complex_type_generator()->element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.description( "Enzymatic mover to glycosylate a pose" )
		.write_complex_type_to_schema( xsd );
}

// Protected methods //////////////////////////////////////////////////////////
void
GlycosyltransferaseMover::perform_reaction(
	core::pose::Pose & input_pose,
	core::uint const site,
	std::string const & cosubstrate )
{
	using namespace core::pose::carbohydrates;

	std::string const atom( get_reactive_site_atom_name( site ) );
	if ( ( atom == "O" ) || ( atom == "-" ) ) {
		// "O" or "-" indicates the O of either Ser or Thr or another unspecified atom name, repsectively.
		// Use the default atom name for the ResidueType being modified.
		glycosylate_pose( input_pose, get_reactive_site_sequence_position( site ), cosubstrate );
	} else {
		glycosylate_pose( input_pose, get_reactive_site_sequence_position( site ), atom, cosubstrate );
	}
}


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

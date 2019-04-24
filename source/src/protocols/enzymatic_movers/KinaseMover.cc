// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/enzymatic_movers/KinaseMover.cc
/// @brief   Method definitions for KinaseMover.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Unit headers
#include <protocols/enzymatic_movers/KinaseMover.hh>
#include <protocols/enzymatic_movers/KinaseMoverCreator.hh>
#include <protocols/enzymatic_movers/EnzymaticMover.hh>

// Project headers
#include <protocols/moves/mover_schemas.hh>

#include <core/chemical/VariantType.hh>
#include <core/pose/variant_util.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.hh>


namespace protocols {
namespace enzymatic_movers {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
KinaseMover::KinaseMover(): EnzymaticMover( "kinases" )
{
	type( "KinaseMover" );
}

// Copy constructor
KinaseMover::KinaseMover( KinaseMover const & object_to_copy ) : EnzymaticMover( object_to_copy )
{}

// Assignment operator
KinaseMover &
KinaseMover::operator=( KinaseMover const & object_to_copy )
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
KinaseMover::register_options()
{
	EnzymaticMover::register_options();
}


// Mover methods
std::string
KinaseMover::get_name() const {
	return mover_name();
}

moves::MoverOP
KinaseMover::clone() const
{
	return utility::pointer::make_shared< KinaseMover >( *this );
}

moves::MoverOP
KinaseMover::fresh_instance() const
{
	return utility::pointer::make_shared< KinaseMover >();
}


void
KinaseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	EnzymaticMover::xml_schema_complex_type_generator()->element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.description( "Enzymatic mover to phosphorylate a pose" )
		.write_complex_type_to_schema( xsd );
}

// Protected methods //////////////////////////////////////////////////////////
void
KinaseMover::perform_reaction(
	core::pose::Pose & input_pose,
	core::uint const site,
	std::string const & /*cosubstrate*/ )
{
	core::pose::add_variant_type_to_pose_residue(
		input_pose,
		core::chemical::PHOSPHORYLATION,
		get_reactive_site_sequence_position( site ) );
}


// Creator methods ////////////////////////////////////////////////////////////
std::string
KinaseMoverCreator::keyname() const {
	return KinaseMover::mover_name();
}

// Return an up-casted owning pointer (MoverOP) to the mover.
protocols::moves::MoverOP
KinaseMoverCreator::create_mover() const {
	return utility::pointer::make_shared< KinaseMover >();
}

void
KinaseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	KinaseMover::provide_xml_schema( xsd );
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that KinaseMover can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, KinaseMover const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace enzymatic_movers
}  // namespace protocols

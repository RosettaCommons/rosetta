// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/enzymatic_movers/DNAMethyltransferaseMover.cc
/// @brief   Method definitions for DNAMethyltransferaseMover.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Unit headers
#include <protocols/enzymatic_movers/DNAMethyltransferaseMover.hh>
#include <protocols/enzymatic_movers/DNAMethyltransferaseMoverCreator.hh>
#include <protocols/enzymatic_movers/EnzymaticMover.hh>

// Project headers
#include <protocols/moves/mover_schemas.hh>

#include <core/chemical/VariantType.hh>
#include <core/pose/variant_util.hh>

// Utility header
#include <utility/tag/XMLSchemaGeneration.hh>

// Basic header
#include <basic/Tracer.hh>

// Construct tracers.
static basic::Tracer TR( "protocols.enzymatic_movers.DNAMethyltransferaseMover" );


namespace protocols {
namespace enzymatic_movers {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
DNAMethyltransferaseMover::DNAMethyltransferaseMover(): EnzymaticMover( "DNA_methyltransferases" )
{
	type( "DNAMethyltransferaseMover" );
}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
DNAMethyltransferaseMover::register_options()
{
	EnzymaticMover::register_options();
}


// Mover methods
std::string
DNAMethyltransferaseMover::get_name() const {
	return mover_name();
}

moves::MoverOP
DNAMethyltransferaseMover::clone() const
{
	return utility::pointer::make_shared< DNAMethyltransferaseMover >( *this );
}

moves::MoverOP
DNAMethyltransferaseMover::fresh_instance() const
{
	return utility::pointer::make_shared< DNAMethyltransferaseMover >();
}


void
DNAMethyltransferaseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	EnzymaticMover::xml_schema_complex_type_generator()->element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.description( "Enzymatic mover to methylate a DNA-containing pose" )
		.write_complex_type_to_schema( xsd );
}

// Protected methods //////////////////////////////////////////////////////////
void
DNAMethyltransferaseMover::perform_reaction(
	core::pose::Pose & input_pose,
	core::uint const site,
	std::string const & /*cosubstrate*/ )
{
	core::uint const seqpos( get_reactive_site_sequence_position( site ) );
	std::string const & name( input_pose.residue( seqpos ).name3() );
	if ( ( name == " DC" ) ) {
		core::pose::add_variant_type_to_pose_residue( input_pose, core::chemical::C5_METHYLATED_NA, seqpos );
	} else {
		TR << "Rosetta can currently only perform 5-methylations of deoxycytidine, ";
		TR << "since other ResidueTypes have not yet been added to the database." << std::endl;
	}
}


// Creator methods ////////////////////////////////////////////////////////////
std::string
DNAMethyltransferaseMoverCreator::keyname() const {
	return DNAMethyltransferaseMover::mover_name();
}

// Return an up-casted owning pointer (MoverOP) to the mover.
protocols::moves::MoverOP
DNAMethyltransferaseMoverCreator::create_mover() const {
	return utility::pointer::make_shared< DNAMethyltransferaseMover >();
}

void
DNAMethyltransferaseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DNAMethyltransferaseMover::provide_xml_schema( xsd );
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that DNAMethyltransferaseMover can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, DNAMethyltransferaseMover const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace enzymatic_movers
}  // namespace protocols

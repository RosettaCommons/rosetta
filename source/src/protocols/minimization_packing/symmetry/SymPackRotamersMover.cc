// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre
/// @author Gutted by Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMoverCreator.hh>

// Project headers

#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.minimization_packing.symmetry.SymPackRotamersMover" );

// Utility Headers
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace minimization_packing {
namespace symmetry {

SymPackRotamersMover::~SymPackRotamersMover() = default;

protocols::moves::MoverOP SymPackRotamersMover::clone() const { return protocols::moves::MoverOP( new  SymPackRotamersMover( *this ) ); }
protocols::moves::MoverOP SymPackRotamersMover::fresh_instance() const { return protocols::moves::MoverOP( new  SymPackRotamersMover ); }

std::string SymPackRotamersMover::get_name() const {
	return mover_name();
}

std::string SymPackRotamersMover::mover_name() {
	return "SymPackRotamersMover";
}

void SymPackRotamersMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_generator_for_pack_rotamers_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Repacks sidechains with user-supplied options, including TaskOperations. \n"
		" NOTE: This Mover is provided for historical support only. The regular PackRotamersMover should handle symmetry transparently and is preferred." )
		.write_complex_type_to_schema( xsd );

	//SymPackRotamersMover description: The symmetric versions of pack rotamers and rotamer trials movers (they take the same tags as asymmetric versions)
}

std::string SymPackRotamersMoverCreator::keyname() const {
	return SymPackRotamersMover::mover_name();
}

protocols::moves::MoverOP
SymPackRotamersMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SymPackRotamersMover );
}

void SymPackRotamersMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymPackRotamersMover::provide_xml_schema( xsd );
}


}
} // moves
} // protocols

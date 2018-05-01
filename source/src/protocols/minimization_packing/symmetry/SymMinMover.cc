// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MinMover.cc
/// @brief
/// @author ashworth
/// @author Gutted by Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/minimization_packing/symmetry/SymMinMoverCreator.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>

// Package headers

#include <basic/Tracer.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.minimization_packing.symmetry.SymMinMover" );

namespace protocols {
namespace minimization_packing {
namespace symmetry {

SymMinMover::~SymMinMover() = default;

protocols::moves::MoverOP SymMinMover::clone() const { return protocols::moves::MoverOP( new  SymMinMover( *this ) ); }
protocols::moves::MoverOP SymMinMover::fresh_instance() const { return protocols::moves::MoverOP( new  SymMinMover ); }

std::string SymMinMover::get_name() const {
	return mover_name();
}

std::string SymMinMover::mover_name() {
	return "SymMinMover";
}

void SymMinMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = MinMover::complex_type_generator_for_min_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Does minimization over sidechain and/or backbone. \n NOTE: This Mover is provided for historical support. The regular MinMover can handle symmetry transparently and should be preferred." )
		.write_complex_type_to_schema( xsd );
	// SymMinMover description: "The symmetric version of min mover (they take the same tags as asymmetric version). Notice that to refine symmetric degrees of freedom, all jumps must be allowed to move with the tag 'jump=ALL'."
}

std::string SymMinMoverCreator::keyname() const {
	return SymMinMover::mover_name();
}

protocols::moves::MoverOP
SymMinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SymMinMover );
}

void SymMinMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymMinMover::provide_xml_schema( xsd );
}



} // symmetry
} // moves
} // protocols

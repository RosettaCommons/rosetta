// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/AddMembraneSpanConstraintCreator.hh
/// @brief      Add membrane span constraint
/// @author     Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_membrane_AddMembraneSpanConstraint_cc
#define INCLUDED_protocols_membrane_AddMembraneSpanConstraint_cc

// Unit Headers
#include <protocols/membrane/AddMembraneSpanConstraint.hh>
#include <protocols/membrane/AddMembraneSpanConstraintCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/util.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/scoring/constraints/MembraneSpanConstraint.hh>

// C++ Headers
#include <cstdlib>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.membrane.AddMembraneSpanConstraint" );

namespace protocols {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
AddMembraneSpanConstraint::AddMembraneSpanConstraint() : protocols::moves::Mover()
{
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
AddMembraneSpanConstraint::AddMembraneSpanConstraint( AddMembraneSpanConstraint const & ) = default;

/// @brief Assignment Operator
AddMembraneSpanConstraint & AddMembraneSpanConstraint::operator = ( AddMembraneSpanConstraint const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new AddMembraneSpanConstraint( *this ) );
}

/// @brief Destructor
AddMembraneSpanConstraint::~AddMembraneSpanConstraint() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
AddMembraneSpanConstraint::clone() const {
	return ( protocols::moves::MoverOP( new AddMembraneSpanConstraint( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
AddMembraneSpanConstraint::fresh_instance() const {
	return protocols::moves::MoverOP( new AddMembraneSpanConstraint() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
AddMembraneSpanConstraint::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
}

/// @brief Create a new copy of this mover
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddMembraneSpanConstraintCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AddMembraneSpanConstraint() );
// XRW TEMP }

/// @brief Return the Name of this mover (as seen by Rscripts)
// XRW TEMP std::string
// XRW TEMP AddMembraneSpanConstraintCreator::keyname() const {
// XRW TEMP  return AddMembraneSpanConstraint::mover_name();
// XRW TEMP }

/// @brief Mover name for Rosetta Scripts
// XRW TEMP std::string
// XRW TEMP AddMembraneSpanConstraint::mover_name() {
// XRW TEMP  return "AddMembraneSpanConstraint";
// XRW TEMP }


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (AddMembraneSpanConstraint)
std::string
AddMembraneSpanConstraint::get_name() const {
	return "AddMembraneSpanConstraint";
}

/// @brief Flip the downstream partner in the membrane
void AddMembraneSpanConstraint::apply( core::pose::Pose & pose ) {
	using namespace core::scoring::constraints;

	ConstraintCOP cst( new core::scoring::constraints::MembraneSpanConstraint( pose ) );
	pose.add_constraint( cst );

}// apply


std::string AddMembraneSpanConstraint::mover_name() {
	return "AddMembraneSpanConstraint";
}

void AddMembraneSpanConstraint::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add membrane span constraint", attlist );
}

std::string AddMembraneSpanConstraintCreator::keyname() const {
	return AddMembraneSpanConstraint::mover_name();
}

protocols::moves::MoverOP
AddMembraneSpanConstraintCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddMembraneSpanConstraint );
}

void AddMembraneSpanConstraintCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddMembraneSpanConstraint::provide_xml_schema( xsd );
}



} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneSpanConstraint_cc

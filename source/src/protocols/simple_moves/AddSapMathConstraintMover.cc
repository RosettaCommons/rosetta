// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddSapMathConstraintMover.cc
/// @brief Mover that adds the SapMathConstraint to the pose
/// @author Brian Coventry (bcov@uw.edu)

// Unit Headers
#include <protocols/simple_moves/AddSapMathConstraintMover.hh>
#include <protocols/simple_moves/AddSapMathConstraintMoverCreator.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/pointer/memory.hh>

// C++ Headers
#include <iostream>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using namespace core;
using namespace core::pack::guidance_scoreterms::sap;

// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {

static basic::Tracer TR( "protocols.simple_moves.AddSapMathConstraintMover" );


// AddSapMathConstraintMover; based on the protocols::moves::Mover basis class
AddSapMathConstraintMover::AddSapMathConstraintMover() :
	protocols::moves::Mover("AddSapMathConstraintMover")
{}


AddSapMathConstraintMover::AddSapMathConstraintMover( AddSapMathConstraintMover const & ot ) :
	protocols::moves::Mover( ot )
{
	*this = ot;
}

AddSapMathConstraintMover &
AddSapMathConstraintMover::operator=( AddSapMathConstraintMover const & ot ) {
	protocols::moves::Mover::operator=( ot );
	cst_ = ot.cst_;
	return *this;
}

void
AddSapMathConstraintMover::apply( core::pose::Pose & pose )
{
	pose.add_constraint( cst_.clone() );
}

moves::MoverOP
AddSapMathConstraintMover::clone() const
{
	return utility::pointer::make_shared< AddSapMathConstraintMover >( *this );
}

moves::MoverOP
AddSapMathConstraintMover::fresh_instance() const
{
	return utility::pointer::make_shared< AddSapMathConstraintMover >( );
}

bool has_option( utility::tag::TagCOP const & tag, std::string const & key ) {
	std::string check = "hasOptionIsBrokenForEmptyStrings";
	return tag->getOption<std::string>( key, check ) != check;
}


void AddSapMathConstraintMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	// Leaving lower_bound blank is valid
	if ( has_option( tag, "lower_bound" ) && !tag->getOption<std::string>( "lower_bound" ).empty() ) {
		set_lower_bound( tag->getOption<Real>( "lower_bound" ) );
	}
	// Leaving upper_bound blank is valid
	if ( has_option( tag, "upper_bound" ) && !tag->getOption<std::string>( "upper_bound" ).empty() ) {
		set_upper_bound( tag->getOption<Real>( "upper_bound" ) );
	}
	if ( tag->hasOption( "penalty_per_unit" ) ) {
		set_penalty_per_unit( tag->getOption<Real>( "penalty_per_unit" ) );
	}

	for ( utility::tag::TagCOP sub_tag : tag->getTags() ) {
		std::string cst_name = sub_tag->getOption<std::string>( "name" );
		Real multiplier = sub_tag->getOption<Real>( "multiplier", Real(1) );

		add_constraint( multiplier, cst_name );
	}
}

void
AddSapMathConstraintMover::add_constraint( Real weight, std::string const & name ) {
	cst_.add_constraint( weight, name );
}

void
AddSapMathConstraintMover::set_upper_bound( Real upper ) {
	cst_.upper_bound( upper );
}

void
AddSapMathConstraintMover::set_lower_bound( Real lower ) {
	cst_.lower_bound( lower );
}

void
AddSapMathConstraintMover::set_penalty_per_unit( Real penalty ) {
	cst_.penalty_per_unit( penalty );
}

std::string AddSapMathConstraintMover::get_name() const {
	return mover_name();
}

std::string AddSapMathConstraintMover::mover_name() {
	return "AddSapMathConstraintMover";
}

void AddSapMathConstraintMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"penalty_per_unit", xsct_real,
		"If the math expression is outside lower_bound or upper_bound, apply this penalty for each unit past the bounds.",
		"1" )
		+ XMLSchemaAttribute::attribute_w_default(
		"lower_bound", xs_string,
		"The lowest value of the math expression that will not incur a penalty. Leaving this blank or not including it will have no lower_bound.",
		"" )
		+ XMLSchemaAttribute::attribute_w_default(
		"upper_bound", xs_string,
		"The highest value of the math expression that will not incur a penalty. Leaving this blank or not including it will have no upper_bound.",
		"" );


	AttributeList subelement_attlist;
	subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "name", xs_string, "Name of constraint added to pose. (Name of AddSapConstraintMover)." )
		+ XMLSchemaAttribute( "multiplier", xsct_real, "Multiply the SapScore by this before adding to the sum." );

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "Add", subelement_attlist, "Specify SapConstraints to be summed." );



	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"A mover that adds the SapMathConstraint to the pose. The SapMathConstraint allows you to create math expressions from previously added"
		" SapConstraints. These SapConstraints are added with the AddSapConstraintMover. The sap_constraint scoreterm must be enabled."
		" Additionally, the SapConstraints this relies on do not need to actually incur a penalty. Setting their penalty_per_sap=0 is what to do"
		" if you only want SapMathConstraints to affect the score."
		" See this paper for more info on sap: Developability index: a rapid in silico tool for the screening of antibody aggregation propensity. "
		" Lauer, et. al. J Pharm Sci 2012",
		attlist, subelements );
}

std::string AddSapMathConstraintMoverCreator::keyname() const {
	return AddSapMathConstraintMover::mover_name();
}

protocols::moves::MoverOP
AddSapMathConstraintMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddSapMathConstraintMover );
}

void AddSapMathConstraintMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddSapMathConstraintMover::provide_xml_schema( xsd );
}


} // moves
} // protocols


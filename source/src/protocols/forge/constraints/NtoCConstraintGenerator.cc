// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/forge/constraints/NtoCConstraintGenerator.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu), Nov 2012

// Unit header
#include <protocols/forge/constraints/NtoCConstraintGenerator.hh>
#include <protocols/forge/constraints/NtoCConstraintGeneratorCreator.hh>

// Package headers
#include <core/pose/Pose.hh>

// Project headers

// utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// boost headers

static basic::Tracer TR( "protocols.forge.constraints.NtoCConstraintGenerator" );

namespace protocols {
namespace forge {
namespace constraints {

// XRW TEMP std::string
// XRW TEMP NtoCConstraintGeneratorCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return NtoCConstraintGenerator::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP NtoCConstraintGeneratorCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new NtoCConstraintGenerator() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP NtoCConstraintGenerator::mover_name()
// XRW TEMP {
// XRW TEMP  return "NtoCConstraintGenerator";
// XRW TEMP }

/// @brief
NtoCConstraintGenerator::NtoCConstraintGenerator():
	RemodelConstraintGenerator(),
	cg_()
{
}

/// @brief
NtoCConstraintGenerator::NtoCConstraintGenerator( Real const dist, Real const coef ):
	RemodelConstraintGenerator(),
	cg_()
{
	cg_.set_weight( coef );
	cg_.set_max_distance( dist );
}

/// @brief
NtoCConstraintGenerator::~NtoCConstraintGenerator() {}

void
NtoCConstraintGenerator::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	cg_.set_weight( tag->getOption< core::Real >( "weight", cg_.weight() ) );
	cg_.set_max_distance( tag->getOption< core::Real >( "dist", cg_.max_distance() ) );
}

// XRW TEMP std::string
// XRW TEMP NtoCConstraintGenerator::get_name() const
// XRW TEMP {
// XRW TEMP  return NtoCConstraintGenerator::mover_name();
// XRW TEMP }

protocols::moves::MoverOP
NtoCConstraintGenerator::fresh_instance() const
{
	return protocols::moves::MoverOP( new NtoCConstraintGenerator() );
}

protocols::moves::MoverOP
NtoCConstraintGenerator::clone() const
{
	return protocols::moves::MoverOP( new NtoCConstraintGenerator( *this ) );
}

/*
/// @brief set weight
void
NtoCConstraintGenerator::set_weight( Real const coef )
{
cg_.set_weight( coef );
}

/// @brief set distance of constraint
void
NtoCConstraintGenerator::set_distance( Real const dist )
{
cg_.set_distance( dist );
}
*/


/// @brief
void
NtoCConstraintGenerator::generate_remodel_constraints( Pose const & pose )
{
	add_constraints( cg_.apply( pose ) );
} //generate constraints

std::string NtoCConstraintGenerator::get_name() const {
	return mover_name();
}

std::string NtoCConstraintGenerator::mover_name() {
	return "NtoCConstraintGenerator";
}

void NtoCConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	RemodelConstraintGenerator::attributes_for_remodel_constraint_generator( attlist );
	//Holds a TerminiConstraintGenerator
	attlist
		+ XMLSchemaAttribute( "weight", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute( "dist", xsct_real, "XRW TO DO" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string NtoCConstraintGeneratorCreator::keyname() const {
	return NtoCConstraintGenerator::mover_name();
}

protocols::moves::MoverOP
NtoCConstraintGeneratorCreator::create_mover() const {
	return protocols::moves::MoverOP( new NtoCConstraintGenerator );
}

void NtoCConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NtoCConstraintGenerator::provide_xml_schema( xsd );
}



} //namespace constraints
} //namespace forge
} //namespace protocols

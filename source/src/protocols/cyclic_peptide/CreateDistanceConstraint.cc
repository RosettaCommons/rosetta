// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
///
/// @file   protocols/cyclic_peptide/CreateDistanceConstraint.cc
/// @brief  Add distance constraints to the current pose conformation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author Yifan Song
///@author modified by Parisa Hosseinzadeh (parisah@uw.edu) to add setters.

#include <protocols/cyclic_peptide/CreateDistanceConstraint.hh>
#include <protocols/cyclic_peptide/CreateDistanceConstraintCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/id/AtomID.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <istream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.CreateDistanceConstraint" );

namespace protocols {
namespace cyclic_peptide {

CreateDistanceConstraint::CreateDistanceConstraint() //:
{}
CreateDistanceConstraint::~CreateDistanceConstraint()= default;

//adding a setter so that people can call the mover within Rosetta commands
void CreateDistanceConstraint::set (
	utility::vector1<Size> const &res1,
	utility::vector1<std::string> const &atom1,
	utility::vector1<Size> const &res2,
	utility::vector1<std::string> const &atom2,
	utility::vector1<std::string> const &cst_function
)
{
	res1_=res1;
	atom1_=atom1;
	res2_=res2;
	atom2_=atom2;
	cst_func_=cst_function;
}

///////////////////////////////////////////////////////////////
/////////////          APPLY FUNCTION         ////////////////
//////////////////////////////////////////////////////////////
//The actual apply function

void CreateDistanceConstraint::apply( core::pose::Pose & pose )
{
	for ( Size i_cst=1; i_cst<=cst_func_.size(); ++i_cst ) {
		if ( cst_func_[i_cst] == "" ) {
			///TODO: use ideal bond length as default parameter

			//Real length, deviation;
			//new core::scoring::func::HarmonicFunc(length, deviation);
			/*
			pose.add_constraint(
			new core::scoring::constraints::AtomPairConstraint( core::id::AtomID(iatom,ires), core::id::AtomID(jatom,jres),
			new core::scoring::func::ScalarWeightedFunc( cst_weight_, new core::scoring::func::Ha( dist, coord_dev_ ) ) ) );
			*/

		} else {
			std::istringstream data(cst_func_[i_cst]);
			std::string func_type;
			data >> func_type;
			core::scoring::func::FuncFactory func_factory;
			core::scoring::func::FuncOP func = func_factory.new_func( func_type );
			func->read_data( data );
			Size atomno1 = pose.residue_type(res1_[i_cst]).atom_index(atom1_[i_cst]);
			Size atomno2 = pose.residue_type(res2_[i_cst]).atom_index(atom2_[i_cst]);
			pose.add_constraint(
				core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint( core::id::AtomID(atomno1,res1_[i_cst]), core::id::AtomID(atomno2,res2_[i_cst]), func ) ) )
			);
		}
	}
	return;
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
CreateDistanceConstraint::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Add" ) {
			res1_.push_back( (*tag_it)->getOption< Size >( "res1" ) );
			atom1_.push_back( (*tag_it)->getOption< std::string >( "atom1" ) );
			res2_.push_back( (*tag_it)->getOption< Size >( "res2" ) );
			atom2_.push_back( (*tag_it)->getOption< std::string >( "atom2" ) );
			cst_func_.push_back( (*tag_it)->getOption< std::string >( "cst_func", "" ) );
		}
	}
	return;
}

moves::MoverOP CreateDistanceConstraint::clone() const { return moves::MoverOP( new CreateDistanceConstraint( *this ) ); }
moves::MoverOP CreateDistanceConstraint::fresh_instance() const { return moves::MoverOP( new CreateDistanceConstraint ); }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CreateDistanceConstraintCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CreateDistanceConstraint );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CreateDistanceConstraintCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return CreateDistanceConstraint::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CreateDistanceConstraint::mover_name()
// XRW TEMP {
// XRW TEMP  return "CreateDistanceConstraint";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CreateDistanceConstraint::get_name() const {
// XRW TEMP  return "CreateDistanceConstraint";
// XRW TEMP }

std::string CreateDistanceConstraint::get_name() const {
	return mover_name();
}

std::string CreateDistanceConstraint::mover_name() {
	return "CreateDistanceConstraint";
}

void CreateDistanceConstraint::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	AttributeList subelement_attlist;
	subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "res1", xsct_non_negative_integer, "First residue to constrain" )
		+ XMLSchemaAttribute::required_attribute( "atom1", xs_string, "First atom to constrain" )
		+ XMLSchemaAttribute::required_attribute( "res2", xsct_non_negative_integer, "Second residue to constrain" )
		+ XMLSchemaAttribute::required_attribute( "atom2", xs_string, "Second atom to constrain" )
		+ XMLSchemaAttribute::required_attribute( "cst_func", xs_string, "Function for distance constraint" );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_simple_subelement( "Add", subelement_attlist, "Specifies a distance constraint between res1,atom1 and res2,atom2" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Adds distance constraints to a pose", attlist, subelements );
}

std::string CreateDistanceConstraintCreator::keyname() const {
	return CreateDistanceConstraint::mover_name();
}

protocols::moves::MoverOP
CreateDistanceConstraintCreator::create_mover() const {
	return protocols::moves::MoverOP( new CreateDistanceConstraint );
}

void CreateDistanceConstraintCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CreateDistanceConstraint::provide_xml_schema( xsd );
}


} // moves
} // protocols

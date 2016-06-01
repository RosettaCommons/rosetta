// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/--class--.hh
/// @brief  --brief-- 
/// @author --name-- (--email--)

// Unit headers
#include <--path--/--class--.hh>
#include --res_sel_creator--
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );


--namespace--

using namespace core::select::residue_selector;

/// @brief Constructor.
///
--class--::--class--():
    ResidueSelector()
{



}

/// @brief Destructor.
///
--class--::~--class--() {}


--class--::--class--( --class-- const & src ):
	TaskOperation(src)

{

}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
--class--::clone() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast<core::select::residue_selector::ResidueSelector>(
			--class--OP( new --class--(*this) )
		)
	);
}

ResidueSelectorOP
--class--Creator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
			--class--OP( new --class-- )
		)
	);
}


std::string --class--::get_name() const
{
	return --class--::class_name();
}

std::string --class--::class_name()
{
	return "--class--";
}

std::string
--class--Creator::keyname() const {
	return --class--::class_name();
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
--class--::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{


}

void --class--::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
    //Syntax Example:
	//using namespace utility::tag;
	//AttributeList attributes;
	//attributes
	//	+ XMLSchemaAttribute::attribute_w_default(  "select_positive_phi",      xs_boolean, "true" )
	//	+ XMLSchemaAttribute::attribute_w_default(  "ignore_unconnected_upper", xs_boolean, "true" );
	//xsd_type_definition_w_attributes( xsd, class_name(), attributes );

}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
--class--Creator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	--class--::provide_xml_schema( xsd );
}



/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResidueSubset
--class--::apply(
	core::pose::Pose const & //pose
) const {

}


--end_namespace--

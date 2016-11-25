// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/PhiSelector.hh
/// @brief  A ResidueSelector that selects alpha-amino acids that are either in the positive phi or negative
/// phi region of Ramachandran space (depending on user preferences).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <core/select/residue_selector/PhiSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.PhiSelector" );


namespace core {
namespace select {
namespace residue_selector {

using namespace core::select::residue_selector;

/// @brief Constructor.
///
PhiSelector::PhiSelector() :
	select_positive_phi_(true),
	ignore_unconnected_upper_(true)
	//TODO -- initialize all vars here.
{}

/// @brief Destructor.
///
PhiSelector::~PhiSelector() {}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
ResidueSelectorOP
PhiSelector::clone() const {
	return ResidueSelectorOP( utility::pointer::dynamic_pointer_cast<ResidueSelector>( PhiSelectorOP( new PhiSelector(*this) ) ) );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResidueSubset
PhiSelector::apply(
	core::pose::Pose const & pose
) const {
	core::Size const nres( pose.size() ); //The number of residues in the pose.

	ResidueSubset selected( nres, false ); //Output, initialized to a vector of "false".

	for ( core::Size i=1; i<=nres; ++i ) { //Loop through all residues.
		if ( !pose.residue_type(i).is_alpha_aa() ) continue; //Skip non-alpha amino acid positions.
		if ( !pose.residue(i).has_lower_connect() || (ignore_unconnected_upper() && !pose.residue(i).has_upper_connect() ) ) continue; //Ignore residues with no lower or upper connection.
		if ( !pose.residue(i).connected_residue_at_resconn( pose.residue_type(i).lower_connect_id() ) ) continue; //Ignore residues with a lower connection, but nothing there.
		if ( ignore_unconnected_upper() && !pose.residue(i).connected_residue_at_resconn( pose.residue_type(i).upper_connect_id() ) ) continue; //Ignore residues with an upper connection, but nothing there.
		if ( pose.phi(i) >= 0 ) { //Select positive phi positions or negative phi positions.
			selected[i] = select_positive_phi();
		} else {
			selected[i] = !select_positive_phi();
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "PhiSelector has selected:" << std::endl;
		for ( core::Size i=1, imax=selected.size(); i<=imax; ++i ) {
			TR.Debug << i << "\t" << "phi:";
			if ( pose.residue(i).type().is_alpha_aa() ) {
				TR.Debug << pose.phi(i);
			} else {
				TR.Debug << "N/A";
			}
			TR.Debug << "\t" << (selected[i] ? "TRUE" : "FALSE") << std::endl;
		}
		TR.Debug.flush();
	}

	return selected;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
PhiSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*datamap*/
) {
	std::string const name( tag->getOption<std::string>( "name" ) );
	if ( TR.visible() ) TR << "Parsing options for PhiSelector \"" << name << "\"." << std::endl;
	set_select_positive_phi( tag->getOption<bool>( "select_positive_phi", select_positive_phi() ) );
	set_ignore_unconnected_upper( tag->getOption<bool>( "ignore_unconnected_upper", ignore_unconnected_upper() ) );
	if ( TR.visible() ) {
		TR << "Set selector to select residues in the "
			<< (select_positive_phi() ? "positive" : "negative") << " phi region of Ramachandran space, "
			<< (ignore_unconnected_upper() ? "ignoring" : "including") << " residues with no upper (C-terminal) connection."
			<< std::endl;
	}
}

/// @brief Get the mover class name.
///
std::string
PhiSelector::get_name() const {
	return PhiSelector::class_name();
}

/// @brief Get the mover class name.
///
std::string
PhiSelector::class_name() {
	return "Phi";
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
PhiSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(
		"select_positive_phi", xsct_rosetta_bool,
		"If true (the default), alpha-amino acids with phi values greater than or "
		"equal to zero are selected. If false, alpha-amino acids "
		"with phi values less than zero are selected.",
		"true" )
		+ XMLSchemaAttribute::attribute_w_default(
		"ignore_unconnected_upper", xsct_rosetta_bool,
		"If true (the default) then C-terminal residues and other residues "
		"with nothing connected at the upper connection are not selected. "
		"If false, then these residues can be selected, depending on their "
		"phi values. Note that anything lacking a lower connection is never selected",
		"true" );

	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The PhiSelector selects alpha-amino acids that are in either the positive "
		"phi or negative phi region of Ramachandran space. Ligands and polymeric "
		"residues that are not alpha-amion acids are never selected. "
		"Alpha-amino acids with no lower connection (or nothing connected at the "
		"lower connection) are also never selected. By default, alpha-amino acids "
		"with no upper connection are not selected, though this can be disabled.",
		attributes );
}

ResidueSelectorOP
PhiSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( utility::pointer::dynamic_pointer_cast<ResidueSelector>( PhiSelectorOP( new PhiSelector ) ) );
}

std::string
PhiSelectorCreator::keyname() const {
	return PhiSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
PhiSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	PhiSelector::provide_xml_schema( xsd );
}


} //core
} //select
} //residue_selector

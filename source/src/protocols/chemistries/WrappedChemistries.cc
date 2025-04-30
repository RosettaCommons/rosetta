// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/WrappedChemistries.cc
/// @brief  XML-able Wrappers for the ChemistryBase-derived chemistries in core
/// @author Rocco Moretti

#include <protocols/chemistries/WrappedChemistries.hh>

#include <protocols/chemistries/util.hh>

#include <core/chemical/AtomRefMapping.hh>

#include <core/chemical/modifications/Reprotonate.hh>

#include <utility/tag/XMLSchemaGeneration.hh>

namespace protocols {
namespace chemistries {

void
WrappedBaseChemistry::apply( core::chemical::MutableResidueType & restype ) {
	runtime_assert( sub_chemistry_ );
	sub_chemistry_->apply( restype );
	set_last_status( sub_chemistry_->get_last_status() );
}

void
WrappedBaseChemistry::apply( core::chemical::MutableResidueType & restype, core::pose::Pose const & ) {
	this->apply( restype );
}

bool
WrappedBaseChemistry::has_additional_output() const {
	runtime_assert( sub_chemistry_ );
	return sub_chemistry_->has_additional_output();
}

core::chemical::MutableResidueTypeOP
WrappedBaseChemistry::get_additional_output() {
	runtime_assert( sub_chemistry_ );
	return sub_chemistry_->get_additional_output();
}

core::chemical::VDVDMapping
WrappedBaseChemistry::get_mapping() const {
	runtime_assert( sub_chemistry_ );
	return sub_chemistry_->get_mapping();
}

////////////////////////////////////////////////////////////////////////////////


ReprotonateChemistry::ReprotonateChemistry() :
	WrappedBaseChemistry( class_name() )
{
	using namespace core::chemical::modifications;
	set_subchemistry( ChemistryBaseOP( new Reprotonate ) );
}

void
ReprotonateChemistry::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &) {
	// Don't need to do anything.
}

std::string
ReprotonateChemistry::class_name() {
	return "ReprotonateChemistry";
}

void
ReprotonateChemistry::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The Reprotonate chemistry will heuristically add and remove hydrogens and charges to change a ResidueType into a form that would be found at neutral pH.",
		attributes );
}

}
}



// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SymmetricalResidueSelector.cc
/// @brief  The SymmetricalResidueSelector selects symmetrical counterparts of the input residues
/// @author Yang Hsia (yhsia@uw.edu)


// Unit headers
#include <core/select/residue_selector/SymmetricalResidueSelector.hh>
#include <core/select/residue_selector/SymmetricalResidueSelectorCreator.hh>

// Protocol Headers

// Core headers
#include <core/select/residue_selector/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/Pose.hh>

// Basic Headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/exit.hh>

// C++ headers

static basic::Tracer TR( "core.select.residue_selector.SymmetricalResidueSelector" );

namespace core {
namespace select {
namespace residue_selector {

SymmetricalResidueSelector::SymmetricalResidueSelector() :
	selector_()
{
}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP SymmetricalResidueSelector::clone() const { return utility::pointer::make_shared< SymmetricalResidueSelector >(*this); }

SymmetricalResidueSelector::SymmetricalResidueSelector(
	core::select::residue_selector::ResidueSelectorCOP const selector ) :
	selector_( selector )
{
}


SymmetricalResidueSelector::~SymmetricalResidueSelector() = default;

ResidueSubset
SymmetricalResidueSelector::apply( core::pose::Pose const & pose ) const
{
	using namespace core::pose::symmetry;
	using namespace core::conformation::symmetry;

	//make sure a selector is set as input
	runtime_assert( selector_ );


	// The asymmetric selection:
	ResidueSubset const subset_non( selector_->apply( pose ) );

	if ( ! core::pose::symmetry::is_symmetric(pose) ) return subset_non;

	// Start the symmetric selection matching the asymmetric:
	ResidueSubset subset_sym( subset_non );

	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);

	// We have to do this in two passes:
	// 1. For each independent position, check if dependents are selected.
	// 2. For each dependent position, check if independent is selected.

	for ( Size passes(1); passes <= 2; ++passes ) {
		for ( Size i = 1, imax(pose.total_residue()); i <= imax; ++i ) {
			if ( TR.Debug.visible() && passes == 1 ) {
				TR.Debug << "Residue " << i << " selected in asymmetric selection: " << (subset_non[i] ? "TRUE" : "FALSE") << std::endl;
			}
			if ( passes == 1 ) {
				// First pass: check the dependent positions and update the independent.
				if ( sym_info->bb_is_independent( i ) ) {
					if ( subset_non[i] ) {
						subset_sym[i] = true;
						continue;
					}
					utility::vector1< core::Size > clones=sym_info->bb_clones( i );
					for ( Size j = 1, jmax(clones.size()); j <= jmax; ++j ) {
						if ( subset_non[ clones[j] ] ) {
							subset_sym[i] = true;
							break;
						}
					}
				}
			} else {
				// Second pass: check the dependent positions of the selected independent positions.
				if ( sym_info->bb_is_independent( i ) ) {
					if ( subset_sym[i] ) {
						utility::vector1< core::Size > clones=sym_info->bb_clones( i );
						for ( Size j = 1, jmax(clones.size()); j <= jmax; ++j ) {
							subset_sym[ clones[j] ] = true;
						}
					}
				}
			}
		}
	}

	//prints all symmetrized residues (output for debugging only)
	if ( TR.Debug.visible() ) {
		TR.Debug << "Symmetrized residues: ";
		for ( Size k = 1; k <= subset_sym.size(); k++ ) {
			if ( subset_sym[k] == 1 ) {
				TR.Debug << k << ",";
			}
		}
		TR.Debug << std::endl;
	}

	return subset_sym;
}

void
SymmetricalResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	//make sure there is only 1 option input
	runtime_assert_string_msg( tag->getTags().size() <= 1, "You can have one or less input tag for SymmetricalResidueSelector!" );

	//if there is a selector, then:
	if ( tag->hasOption( "selector" ) ) {
		std::string const selectorname = tag->getOption< std::string >( "selector" );
		//make sure it is a valid residue selector
		try {
			set_selector( core::select::residue_selector::get_residue_selector(selectorname, data) );
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selectorname << "' from the Datamap from SymmetricalResidueSelector.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  error_msg.str() );
		}
		debug_assert( selector_ );
		TR << "Using residue selector " << selectorname << std::endl;
	} else {
		/*
		//Shouldn't need this since there can only be 1 option tag?
		else if ( tag->getTags().size() ) {
		for ( utility::vector0< utility::tag::TagCOP >::const_iterator t = tag->getTags().begin();
		t != tag->getTags().end(); ++t ) {
		ResidueSelectorCOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
		(*t)->getName(),
		(*t),
		data );
		if ( rs ) {
		debug_assert( !selector_ );
		set_selector( rs );
		}
		}
		} */
		//if there is not a selector tag option, then:
		std::stringstream ss;
		ss << "SymmetricalResidueSelector requires a selector to be set either through the \"selector\" option." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  ss.str() );
	}
	debug_assert( selector_ );
}

void
SymmetricalResidueSelector::set_selector(
	core::select::residue_selector::ResidueSelectorCOP const selector )
{
	selector_ = selector;
}

std::string SymmetricalResidueSelector::get_name() const {
	return SymmetricalResidueSelector::class_name();
}

std::string SymmetricalResidueSelector::class_name() {
	return "SymmetricalResidue";
}

void SymmetricalResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes_for_parse_residue_selector(attributes, "selector");

	xsd_type_definition_w_attributes_and_optional_subselector(
		xsd, class_name(),
		"The SymmetricalResidueSelector, when given a selector, will return all symmetrical "
		"copies (including the original selection) of those residues. While the packer is "
		"symmetry aware, not all filters are. This selector is useful when you need to "
		"explicitly give residue numbers but you are not sure which symmetry subunit you need.",
		attributes );
}


ResidueSelectorOP
SymmetricalResidueSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared< SymmetricalResidueSelector >();
}

std::string
SymmetricalResidueSelectorCreator::keyname() const {
	return SymmetricalResidueSelector::class_name();
}

void
SymmetricalResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymmetricalResidueSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core

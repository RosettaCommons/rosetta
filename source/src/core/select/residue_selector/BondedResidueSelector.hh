// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/BondedResidueSelector.hh
/// @brief  Selects residues with chemical bonds to the passed input residues
/// @author Sharon Guffy (guffy@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_BondedResidueSelector_HH
#define INCLUDED_core_select_residue_selector_BondedResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/BondedResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility Headers
//#include <utility/tag/Tag.fwd.hh>
//#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// C++ headers
#include <set>

namespace core {
namespace select {
namespace residue_selector {

/// @brief The BondedResidueSelector selects residues that are connected to some
/// set of residues (input_set) by a chemical bond. The input_set can be directly
/// set using a set of residue positions or by a ResidueSelector.
///
/// Residues in the input_set are included in the final selection.
///
///
class BondedResidueSelector : public ResidueSelector {
public:
	// derived from base class
	BondedResidueSelector();
	BondedResidueSelector( std::set<core::Size> const & input_set );
	BondedResidueSelector ( BondedResidueSelector const & other );
	virtual ~BondedResidueSelector();

	virtual ResidueSelectorOP clone() const;
	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();


	//unit-specific
	/**
	* @brief adds a ResidueSelector
	*/
	void set_input_set( std::set<core::Size> const & input_set );
	void set_input_set( std::string const & input_set_str );
	void set_input_set_selector(  core::select::residue_selector::ResidueSelectorCOP rs );

	//Getters (for copy constructor)
	std::set< core::Size > input_set() const;
	core::select::residue_selector::ResidueSelectorCOP input_set_selector() const;
	std::string input_set_string() const;
	bool input_set_defined() const;
	bool input_set_selector_defined() const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: //methods
	void get_input_set( core::pose::Pose const &, ResidueSubset &, std::set< core::Size > &) const;





private: // data members
	// data in focus and focus_string will be stitched together.
	// think of either-or behavior also between set and string
	std::set< core::Size > input_set_;
	std::string input_set_str_;

	// focus residues may be selected directly be another ResidueSelector
	core::select::residue_selector::ResidueSelectorCOP input_set_selector_;

	// have we defined an input set
	bool input_set_defined_;

	// is the input set defined by a selector?
	bool use_input_set_selector_;
};

} //namespace residue_selector
} //namespace select
} //namespace core

#endif

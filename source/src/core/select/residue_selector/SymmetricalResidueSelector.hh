// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/SymmetricalResidueSelector.hh
/// @brief  The SymmetricalResidueSelector selects symmetrically generated residues of a given ResidueSelector input
/// @author Yang Hsia (yhsia@uw.edu)

#ifndef INCLUDED_core_select_residue_selector_SymmetricalResidueSelector_HH
#define INCLUDED_core_select_residue_selector_SymmetricalResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/SymmetricalResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace core {
namespace select {
namespace residue_selector {

/// @brief The SymmetricalResidueSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which are symmetrically generated residues of the input list.
class SymmetricalResidueSelector : public ResidueSelector {
public:
	SymmetricalResidueSelector();

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	SymmetricalResidueSelector(
		core::select::residue_selector::ResidueSelectorCOP const selector );

	virtual ~SymmetricalResidueSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );

	virtual
	std::string
	get_name() const;

	static std::string class_name();

	//unit-specific
	void set_selector( core::select::residue_selector::ResidueSelectorCOP const selector );

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // data members
	core::select::residue_selector::ResidueSelectorCOP selector_;
};

} //namespace residue_selector
} //namespace select
} //namespace core


#endif

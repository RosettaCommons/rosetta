// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/residue_selectors/LigandMetalContactSelectorCreator.hh
/// @brief This residue selector takes a selector or residue number of a ligand and returns any residues in contact with metal atoms in the ligand.
/// @author Allison Watwood (acw538@msstate.edu)

#ifndef INCLUDED_protocols_residue_selectors_LigandMetalContactSelectorCreator_HH
#define INCLUDED_protocols_residue_selectors_LigandMetalContactSelectorCreator_HH

// Package headers
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

namespace protocols {
namespace residue_selectors {

class LigandMetalContactSelectorCreator : public  core::select::residue_selector::ResidueSelectorCreator {
public:
	core::select::residue_selector::ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

} //protocols
} //residue_selectors

#endif //INCLUDED_protocols_residue_selectors_LigandMetalContactSelectorCreator_HH

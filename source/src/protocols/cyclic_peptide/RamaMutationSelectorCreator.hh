// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/RamaMutationSelectorCreator.hh
/// @brief Selects positions that would have a rama_prepro score below a given threshold IF mutated to a given residue type.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_RamaMutationSelectorCreator_HH
#define INCLUDED_protocols_cyclic_peptide_RamaMutationSelectorCreator_HH

// Package headers
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

namespace protocols {
namespace cyclic_peptide {

class RamaMutationSelectorCreator : public  core::select::residue_selector::ResidueSelectorCreator {
public:
	core::select::residue_selector::ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_RamaMutationSelectorCreator_HH

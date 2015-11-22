// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ClashBasedRepackShellSelectorCreator.hh
/// @brief  The ClashBasedRepackShellSelector identifies all residues that clash with at least one rotamer of a design position
/// @details Since this ResidueSelector is located in a different namespace, it needs a separate ResidueSelectorCreator.
/// @author Noah Ollikainen (nollikai@gmail.com)
/// @author Roland A. Pache, PhD
/// @author Vikram K. Mulligan, PhD (vmullig@uw.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_ClashBasedRepackShellSelectorCreator_HH
#define INCLUDED_core_pack_task_residue_selector_ClashBasedRepackShellSelectorCreator_HH

// Unit headers
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class ClashBasedRepackShellSelectorCreator : public core::select::residue_selector::ResidueSelectorCreator {
public:
	virtual core::select::residue_selector::ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
};

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif

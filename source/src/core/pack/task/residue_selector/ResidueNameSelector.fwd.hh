// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueNameSelector.fwd.hh
/// @brief  Forward declaration of a class that selects residues using their names
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_ResidueNameSelector_FWD_HH
#define INCLUDED_core_pack_task_residue_selector_ResidueNameSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.fwd.hh>


namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class ResidueNameSelector;

typedef utility::pointer::shared_ptr< ResidueNameSelector > ResidueNameSelectorOP;
typedef utility::pointer::shared_ptr< ResidueNameSelector const > ResidueNameSelectorCOP;

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif

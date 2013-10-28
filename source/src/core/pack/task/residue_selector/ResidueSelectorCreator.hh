// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueSelectorCreator.hh
/// @brief  Class for instantiating a particular ResidueSelector
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_ResidueSelectorCreator_HH
#define INCLUDED_core_pack_task_residue_selector_ResidueSelectorCreator_HH

// Package headers
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class ResidueSelectorCreator : public utility::pointer::ReferenceCount {
public:
	virtual ResidueSelectorOP create_residue_selector() const = 0;
	virtual std::string keyname() const = 0;
};


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif

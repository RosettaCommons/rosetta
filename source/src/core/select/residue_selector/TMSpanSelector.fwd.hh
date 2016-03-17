// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/select/residue_selector/TMSpanSelector.fwd.hh
/// @brief Select residues within given transmembrane spans in a membrane protein
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_TMSpanSelector_fwd_hh
#define INCLUDED_core_select_residue_selector_TMSpanSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace select {
namespace residue_selector {

class TMSpanSelector;

typedef utility::pointer::shared_ptr< TMSpanSelector > TMSpanSelectorOP;
typedef utility::pointer::shared_ptr< TMSpanSelector const > TMSpanSelectorCOP;

} // residue_selector
} // select
} // core

#endif // INCLUDED_core_select_residue_selector_TMSpanSelector_fwd_hh


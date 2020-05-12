// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pdbinfo_manipulations/AlignPDBInfoToSequences.fwd.hh
/// @brief Realign poses to sequences after losing pdb_info
/// @author Dan Farrell (danpf@uw.edu)

#ifndef INCLUDED_protocols_pdbinfo_manipulations_AlignPDBInfoToSequences_fwd_hh
#define INCLUDED_protocols_pdbinfo_manipulations_AlignPDBInfoToSequences_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace pdbinfo_manipulations {

class AlignPDBInfoToSequences;

using AlignPDBInfoToSequencesOP = utility::pointer::shared_ptr< AlignPDBInfoToSequences >;
using AlignPDBInfoToSequencesCOP = utility::pointer::shared_ptr< AlignPDBInfoToSequences const >;

} //pdbinfo_manipulations
} //protocols

#endif //INCLUDED_protocols_pdbinfo_manipulations_AlignPDBInfoToSequences_fwd_hh

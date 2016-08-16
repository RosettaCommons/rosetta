// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter.fwd.hh
/// @brief Forward declarations and owning pointers for the OversaturatedHbondAcceptorFilter.
/// @details This filter flags poses containing more than two hydrogen bonds to an oxygen atom, a common pathology that results from Rosetta's
/// pairwise-decomposible scorefunction, which can't penalize excessive hydrogen bonds.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_protocols_cyclic_peptide_OversaturatedHbondAcceptorFilter_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_OversaturatedHbondAcceptorFilter_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace cyclic_peptide {

class OversaturatedHbondAcceptorFilter;

typedef utility::pointer::shared_ptr< OversaturatedHbondAcceptorFilter > OversaturatedHbondAcceptorFilterOP;
typedef utility::pointer::shared_ptr< OversaturatedHbondAcceptorFilter const > OversaturatedHbondAcceptorFilterCOP;

} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_OversaturatedHbondAcceptorFilter_fwd_hh

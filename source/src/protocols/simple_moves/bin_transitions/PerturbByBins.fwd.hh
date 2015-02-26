// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/bin_transitions/PerturbByBins.cc
/// @brief  Fwd declarations for the PerturbByBins mover.  This mover takes a stretch of backbone and perturbs its mainchain torsions based
/// on the probabilities of transitions from one torsion bin to another.
/// @details Bin transitions are read from database files.  The algorithm is: set the first residue based on the probability of a residue
/// being in a bin.  Set subsequent residues based on the probability of a residue being in a bin given that the previous residue is in
/// a particular bin.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_bin_transitions_PerturbByBins_fwd_hh
#define INCLUDED_protocols_simple_moves_bin_transitions_PerturbByBins_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

// Package headers

namespace protocols {
	namespace simple_moves {
		namespace bin_transitions {

			class PerturbByBins;
			typedef utility::pointer::shared_ptr< PerturbByBins > PerturbByBinsOP;
			typedef utility::pointer::shared_ptr< PerturbByBins const > PerturbByBinsCOP;
			typedef utility::vector1<PerturbByBinsOP> PerturbByBinsOPs;
			typedef utility::vector1<PerturbByBinsCOP> PerturbByBinsCOPs;

		} // bin_transitions
	} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_bin_transitions_PerturbByBins_fwd_hh

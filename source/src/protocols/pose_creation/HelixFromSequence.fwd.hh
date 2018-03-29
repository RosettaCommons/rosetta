// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/HelixFromSequence.fwd.hh
/// @brief Generates a (transmembrane) helix from a fasta file containing the sequence
/// @author Julia Koehler Leman (julia.koehler1982@gmail.com)


#ifndef INCLUDED_protocols_membrane_HelixFromSequence_fwd_hh
#define INCLUDED_protocols_membrane_HelixFromSequence_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace membrane {

class HelixFromSequence;

typedef utility::pointer::shared_ptr< HelixFromSequence > HelixFromSequenceOP;
typedef utility::pointer::shared_ptr< HelixFromSequence const > HelixFromSequenceCOP;



} //protocols
} //membrane


#endif //INCLUDED_protocols_membrane_HelixFromSequence_fwd_hh






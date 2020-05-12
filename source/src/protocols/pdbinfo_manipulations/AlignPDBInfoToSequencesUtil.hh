// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesUtil.hh
/// @brief Helper functions for AlignPDBInfoToSequences
/// @author Dan Farrell (danpf@uw.edu)

#ifndef INCLUDED_protocols_pdbinfo_manipulations_AlignPDBInfoToSequencesUtil_HH
#define INCLUDED_protocols_pdbinfo_manipulations_AlignPDBInfoToSequencesUtil_HH

// Core headers
#include <core/sequence/Sequence.hh>

// Standard library headers
#include <string>

namespace protocols {
namespace pdbinfo_manipulations {

///@brief Take two sequences that have been aligned, and their parent sequences
///       and align/pad them with '-' as shown below:
///
///       given to align from SWAlign:
///       query: CLGIGSCNDFAGCGYAVVCFW
///       target: XXXXXXXXXXXXXXCLGIGSCNDFAGCGXXXXXXYAVVCFWXXXXX
///
///       from: (aka a result from SWAlign)
///       query   1 CLGIGSCNDFAGCG------YAVVCFW
///       target 15 CLGIGSCNDFAGCGXXXXXXYAVVCFW
///
///       to:
///       query   1 --------------CLGIGSCNDFAGCG------YAVVCFW-----
///       target  1 XXXXXXXXXXXXXXCLGIGSCNDFAGCGXXXXXXYAVVCFWXXXXX
void
pad_sequences(
	std::string const & full_seq_1,
	core::sequence::SequenceOP const & aln_seq_1,
	std::string const & full_seq_2,
	core::sequence::SequenceOP const & aln_seq_2);

}
}

#endif //protocols_pdbinfo_manipulations_AlignPDBInfoToSequencesUtil_HH

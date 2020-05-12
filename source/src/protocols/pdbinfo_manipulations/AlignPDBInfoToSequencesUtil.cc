// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesUtil.cc
/// @brief helper functions for AlignPDBInfoToSequences
/// @author Dan Farrell (danpf@uw.edu)

// Core headers
#include <protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesUtil.hh>
#include <core/sequence/Sequence.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/string_util.hh>

// Standard library headers
#include <string>


static basic::Tracer TR( "protocols.pdbinfo_manipulations.AlignPDBInfoToSequencesUtil" );

namespace protocols {
namespace pdbinfo_manipulations {

void
pad_sequences(
	std::string const & full_seq_1,
	core::sequence::SequenceOP const & aln_seq_1,
	std::string const & full_seq_2,
	core::sequence::SequenceOP const & aln_seq_2)
{
	core::Size const max_start = std::max(aln_seq_1->start(), aln_seq_2->start());

	std::string const gapless_aln_1(utility::remove_from_string(aln_seq_1->sequence(), "-"));
	std::string const gapless_aln_2(utility::remove_from_string(aln_seq_2->sequence(), "-"));

	core::Size const missing_seq_1_size(full_seq_1.size()-gapless_aln_1.size());
	core::Size const missing_seq_2_size(full_seq_2.size()-gapless_aln_2.size());

	std::string const seq_1_final_pad(
		missing_seq_1_size > 0 ? full_seq_1.substr( gapless_aln_1.size() + (aln_seq_1->start()-1), missing_seq_1_size) : "");

	std::string const seq_2_final_pad(missing_seq_2_size > 0 ? full_seq_2.substr( gapless_aln_2.size() + (aln_seq_2->start()-1), missing_seq_2_size) : "");

	std::string const aligned_seq_1(
		std::string(max_start-aln_seq_1->start(), '-')
		+ full_seq_1.substr(0, aln_seq_1->start()-1)
		+ std::string(aln_seq_1->sequence()).erase( aln_seq_1->sequence().find_last_not_of("-") + 1)
		+ seq_1_final_pad
	);

	std::string const aligned_seq_2(
		std::string(max_start-aln_seq_2->start(), '-')
		+ full_seq_2.substr(0, aln_seq_2->start()-1)
		+ std::string(aln_seq_2->sequence()).erase( aln_seq_2->sequence().find_last_not_of("-") + 1)
		+ seq_2_final_pad
	);

	core::Size const max_final_size(std::max(aligned_seq_1.size(), aligned_seq_2.size()));

	aln_seq_1->sequence(aligned_seq_1 + std::string(max_final_size-aligned_seq_1.size(), '-'));
	aln_seq_2->sequence(aligned_seq_2 + std::string(max_final_size-aligned_seq_2.size(), '-'));
	aln_seq_1->start(1);
	aln_seq_2->start(1);
}


}
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/sequence/sequtil.hh
/// @brief small bundle of utilities for dealing with sequences
/// @author James Thompson
/// @author Sergey Lyskov

#ifndef INCLUDED_core_sequence_util_hh
#define INCLUDED_core_sequence_util_hh

// C/C++ headers

// Utility headers
#include <utility/vector1.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/SequenceMapping.fwd.hh>

// Package headers
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <core/sequence/DerivedSequenceMapping.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace sequence {

extern utility::vector1< char > spacers;

/// @brief Populates the non-null vector <alignments> with all of the sequence
/// alignments found in <files>. Each alignment is required to have
/// format <format>.
void read_all_alignments(const std::string& format,
	const utility::vector1<std::string>& files,
	utility::vector1<SequenceAlignment>* alignments);

/// @brief helper function for reading a SequenceMapping from an alignment
/// file.
void
read_alignment_file(
	std::string const & filename,
	std::string & seq1,
	std::string & seq2,
	core::sequence::DerivedSequenceMapping & mapping // from numbering in sequence 1 to numbering in sequence 2
);

/// @brief Read in sequences from a fasta-formatted file.
utility::vector1< SequenceOP > read_fasta_file( std::string const & filename );
utility::vector1< std::string > read_fasta_file_str( std::string const & filename );
std::string read_fasta_file_return_str( std::string const & filename );

///  @brief read sequence from particular section of fasta file (comment starting with '> section'), terminate with failure if section not found
///         Note: section detection string is case insensitive
std::string read_fasta_file_section(std::string const & filename, std::string const &section);

/// @brief Return a string of concatenated SequenceCOP sequences
std::string get_concatenated_sequence( utility::vector1< SequenceCOP > const & fasta_sequences );

/// @brief Read fasta file and concatenate sequences
std::string read_fasta_file_and_concatenate( std::string const & filename );

/// @brief Read in a SequenceMapping from a file. File
/// format is super-simple, it just contains single
/// lines like this that claim that residue resi and
/// resj are aligned: resi resj
core::sequence::DerivedSequenceMapping simple_mapping_from_file( std::string const & filename );

utility::vector1< SequenceAlignment > read_aln(
	std::string const & format,
	std::string const & filename
);

utility::vector1< SequenceOP > seqs_from_cmd_lines();

/// @brief read generalized alignment format.
utility::vector1< SequenceAlignment > read_general_aln(
	std::istream & input
);

utility::vector1< SequenceAlignment > read_grishin_aln_file(
	std::string const & filename
);

utility::vector1< SequenceAlignment > read_general_aln_file(
	std::string const & filename
);

// @brief returns the number of correctly aligned positions in candidate_aln
// relative to true_aln.
core::Size n_correctly_aligned_positions(
	SequenceAlignment & candidate_aln,
	SequenceAlignment & true_aln
);

/// @brief takes the sequences in the provided vector1 and makes them match
/// the alignment in aln_to_steal by matching gaps. This assumes that the
/// ungapped sequences at index j in the vector1< SequenceOP > match the
/// ungapped sequences at index j in aln_to_steal.
SequenceAlignment steal_alignment(
	SequenceAlignment aln_to_steal,
	utility::vector1< SequenceOP > seqs
);

/// @brief Constructs a SequenceAlignment from the given SequenceMapping and
/// the two sequences.
SequenceAlignment mapping_to_alignment(
	core::id::SequenceMapping const & mapping,
	SequenceOP seq1,
	SequenceOP seq2
);

/// @brief Assuming that map1 maps sequence A to sequence B, and map2 maps
/// sequence B to sequence C, this function returns the SequenceMapping
/// representing the direct map of sequence A to sequence C.
core::id::SequenceMapping transitive_map(
	core::id::SequenceMapping const & map1,
	core::id::SequenceMapping const & map2
);

/// @brief Generates a mapping of sequence 1 onto sequence 2 using dynamic
/// programming with a simple scoring framework.
core::id::SequenceMapping map_seq1_seq2(
	core::sequence::SequenceOP seq1,
	core::sequence::SequenceOP seq2
);

/// @brief Generate a naive sequence alignment between two sequences.
core::sequence::SequenceAlignment align_naive(
	core::sequence::SequenceOP seq1,
	core::sequence::SequenceOP seq2
);

core::sequence::SequenceAlignment align_poses_naive(
	core::pose::Pose & pose1,
	core::pose::Pose & pose2
);

utility::vector1< Real >
get_maximum_scores(
	core::sequence::ScoringSchemeOP ss,
	core::sequence::SequenceOP seq
);

core::sequence::SequenceAlignment
alignment_from_pose(
	core::pose::Pose & pose
);

void alignment_into_pose(
	core::sequence::SequenceAlignment const & aln,
	core::pose::Pose & pose
);
/*
core::Real
calpha_superimpose_via_alignment(
core::pose::Pose & mod_pose,
core::pose::Pose const & ref_pose,
core::sequence::SequenceAlignment const & aln
); */

core::Real
calpha_superimpose_with_mapping(
	core::pose::Pose & mod_pose,
	core::pose::Pose const & ref_pose,
	core::id::SequenceMapping const & mapping // mod_pose -> ref_pose
);

utility::vector1< core::Size >
strip_spacers( std::string & sequence );

} // sequence
} // core

#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/enzymes/consensus_sequence_parsers.hh
/// @brief   Declarations for helper functions for the parsing of enzyme consensus sequences.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_enzymes_consensus_sequence_parsers_HH
#define INCLUDED_core_enzymes_consensus_sequence_parsers_HH

// Utility header
#include <utility/vector1.hh>

// C++ header
#include <string>


namespace core {
namespace enzymes {

/// @brief  Parse a peptide consensus sequence and return a list of AA residue 3-letter codes.
utility::vector1< utility::vector1< std::string > > get_3_letter_codes_from_peptide_consensus_sequence(
		std::string const & sequence );

/// @brief  Parse a nucleic acid consensus sequence and return a list of NA residue codes.
utility::vector1< utility::vector1< std::string > > get_codes_from_NA_consensus_sequence(
		std::string const & sequence );

/// @brief  Parse a saccharide consensus sequence and return a list of monosaccharide residue codes.
utility::vector1< utility::vector1< std::string > > get_codes_from_saccharide_consensus_sequence(
		std::string const & sequence );

}  // namespace enzymes
}  // namespace core

#endif  // INCLUDED_core_enzymes_consensus_sequence_parsers_HH

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/enzymes/consensus_sequence_parsers.cc
/// @brief   Definitions for helper functions for the parsing of enzyme consensus sequences.
/// @author  Labonte <JWLabonte@jhu.edu>

// Unit header
#include <core/enzymes/consensus_sequence_parsers.hh>

// Project header
#include <core/types.hh>
#include <core/chemical/AA.hh>


namespace core {
namespace enzymes {

// Parse a peptide consensus sequence and return a list of AA residue 3-letter codes.
/// @note  A vector of vectors is returned, because each position may have options and not a specific residue.
utility::vector1< utility::vector1< std::string > >
get_3_letter_codes_from_peptide_consensus_sequence( std::string const & sequence )
{
	using namespace std;
	using namespace utility;
	using namespace core::chemical;

	vector1< vector1< string > > residues_in_consensus;  // a vector of vectors because a position may have options

	Size const n_char( sequence.size() );
	bool in_parentheses( false );
	vector1< chemical::AA > AA_options;
	for ( uint char_pos( 0 ); char_pos < n_char; ++char_pos ) {
		char const one_letter_code( sequence[ char_pos ] );
		if ( one_letter_code == '(' ) {
			in_parentheses = true;
			AA_options.clear();
			continue;
		} else if ( one_letter_code == ')' ) {
			in_parentheses = false;

			// We have finished collecting a list of optional AA residues.
			// Convert them into 3-letter codes and store them as options for the current position in the
			// consensus sequence.
			Size const n_options( AA_options.size() );
			vector1< string > three_letter_code_options( n_options );
			for ( uint j( 1 ); j <= n_options; ++j ) {
				three_letter_code_options[ j ] = name_from_aa( AA_options[ j ] );
			}
			residues_in_consensus.push_back( three_letter_code_options );
			continue;
		} else if ( one_letter_code == '/' ) {
			continue;
		}
		if ( in_parentheses ) {
			// If we are within a set of parentheses, every 1-letter code represents an optional AA residue.
			AA_options.push_back( aa_from_oneletter_code( one_letter_code ) );
		} else {
			if ( one_letter_code == 'X' ) {  // TODO: Add other 1-letter codes specifying multiple options.
				Size const n_options( 20 );  // There are 20 canonical AAs in the enum.
				vector1< string > three_letter_code_options( n_options );
				for ( uint j( 1 ); j <= n_options; ++j ) {
					three_letter_code_options[ j ] = name_from_aa( static_cast< chemical::AA >( j ) );
				}
				residues_in_consensus.push_back( three_letter_code_options );
			} else {
				// There is only one option for this position.  Convert to a 3-letter code and push back.
				residues_in_consensus.push_back(
					vector1< string >( 1, name_from_aa( aa_from_oneletter_code( one_letter_code ) ) ) );
			}
		}
	}
	return residues_in_consensus;
}

// Parse a nucleic acid consensus sequence and return a list of NA residue codes.
/// @note  A vector of vectors is returned, because each position may have options and not a specific residue.
utility::vector1< utility::vector1< std::string > >
get_codes_from_NA_consensus_sequence( std::string const & /*sequence*/ )
{
	using namespace std;
	using namespace utility;
	return vector1< vector1< string> >( 1, vector1< string >( 1, "XXX" ) );  // TODO: Fill in.
}

// Parse a saccharide consensus sequence and return a list of monosaccharide residue codes.
/// @note  A vector of vectors is returned, because each position may have options and not a specific residue.
utility::vector1< utility::vector1< std::string > >
get_codes_from_saccharide_consensus_sequence( std::string const & /*sequence*/ )
{
	using namespace std;
	using namespace utility;
	return vector1< vector1< string> >( 1, vector1< string >( 1, "XXX" ) );  // TODO: Fill in.
}

}  // namespace enzymes
}  // namespace core

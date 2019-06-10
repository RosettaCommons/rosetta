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

// This helper function for get_3_letter_codes_from_peptide_consensus_sequence() returns a list of AA residue 3-letter
// codes or a 1-element vector with a single-character "punctuation" symbol of "(", "/", or ")".
// It is intentionally not forward-declared.
// <char_pos> is passed as a non-const reference, because the one-letter code "X", if followed by square brackets,
// counts as a single character with the bracket contents specifying its identity, and the next character parsed by
// get_3_letter_codes_from_peptide_consensus_sequence() needs to be the one following the "]".
utility::vector1< std::string >
get_AA_3_letter_codes_or_punctuation_from_1_letter_code_at_position( std::string const & sequence, uint & char_pos )
{
	using namespace std;
	using namespace utility;
	using namespace core::chemical;

	char const first_char( sequence[ char_pos ] );
	if ( ( first_char == '(' ) || ( first_char == '/' ) || ( first_char == ')' ) ) {
		return vector1< string >( { string( 1, first_char ) } );
	} else if ( first_char == 'B' ) {
		return vector1< string >( { "ASP", "ASN" } );
	} else if ( first_char == 'Z' ) {
		return vector1< string >( { "GLU", "GLN" } );
	} else if ( first_char == 'J' ) {
		return vector1< string >( { "LEU", "ILE" } );
	} else if ( first_char == 'O' ) {
		return vector1< string >( { "PYL" } );
	} else if ( first_char == 'U' ) {
		return vector1< string >( { "SEC" } );
	} else if ( first_char == 'X' ) {
		char next_char( sequence[ char_pos + 1 ] );
		if ( next_char == '[' ) {  // We are signaling for a specific NCAA.
			++char_pos;
			bool in_brackets( true );
			string specific_name;
			while ( in_brackets ) {
				next_char = sequence[ char_pos + 1 ];
				if ( next_char == ']' ) {
					in_brackets = false;
				} else {
					specific_name += next_char;
				}
				++char_pos;
			}
			return vector1< string >( { specific_name } );
		} else {  // We are signaling for any of the CAAs.
			Size const n_CAAs( 20 );  // There are 20 canonical AAs in the enum.
			vector1< string > all_CAAs( n_CAAs );
			for ( uint i( 1 ); i <= n_CAAs; ++i ) {
				all_CAAs[ i ] = name_from_aa( static_cast< chemical::AA >( i ) );
			}
			return all_CAAs;
		}
	} else {
		return vector1< string >( { name_from_aa( aa_from_oneletter_code( first_char ) ) } );
	}

}

// Parse a peptide consensus sequence and return a list of AA residue 3-letter codes.
/// @details  This parser recognizes the IUPAC-approved one-letter codes B, J, O, U, and Z, which code for Asx, Xle,
/// Pyl, Sec, and Glx, respectively.  X alone is recognized to be any of the 20 canonical amino acids; X followed by
/// square brackets specifies a single non-canonical amino acid by 3-letter code. For example, X[SEP] specifies
/// phosphoserine.  Parentheses are used to specify multiple possible residue types at that site, separated by forward
/// slashes, e.g., A/G) specifies either Ala or Gly at that position.
/// @note     A vector of vectors is returned, because each position may have options and not a specific residue.
utility::vector1< utility::vector1< std::string > >
get_3_letter_codes_from_peptide_consensus_sequence( std::string const & sequence )
{
	using namespace std;
	using namespace utility;

	vector1< vector1< string > > residues_in_consensus;  // a vector of vectors because a position may have options

	Size const n_char( sequence.size() );
	bool in_parentheses( false );
	vector1< string > AA_options;
	for ( uint char_pos( 0 ); char_pos < n_char; ++char_pos ) {
		vector1< string > const codes_or_punctuation(
			get_AA_3_letter_codes_or_punctuation_from_1_letter_code_at_position( sequence, char_pos ) );
		char const first_char( codes_or_punctuation[ 1 ][ 0 ] );
		if ( first_char == '(' ) {
			in_parentheses = true;
		} else if ( first_char == '/' ) {
			continue;
		} else if ( first_char == ')' ) {
			in_parentheses = false;

			// We have finished collecting a list of optional AA residues.
			// Store them as options for the current position in the consensus sequence.
			residues_in_consensus.push_back( AA_options );
			AA_options.clear();
		} else {
			if ( in_parentheses ) {
				// If we are within a set of parentheses, every 1-letter code represents (an) optional AA residue(s).
				AA_options.append( codes_or_punctuation );
			} else {
				residues_in_consensus.push_back( codes_or_punctuation );
			}
		}
	}
	return residues_in_consensus;
}

// TODO: copied from Slack convo. w/ Andy:
// B = ‘not A’; D = ‘not C’; H = ‘not G’; V = ‘not U’
// R = AG, Y = CU
// S = GC, W = AU (weak/strong)
// K = GU, M = AC (idk why)
// N = ACGU

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

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MatrixScoringScheme.cc
/// @brief class definition for a given scoring scheme for an alignment.
/// @details Simply based on comparing single characters from two protein
/// sequences, along with affine gap penalties of the form penalty = A + Bk, where
/// A represents the penalty for starting a gap, and B represents the penalty for
/// extending a previously opened gap by k characters.
/// @author James Thompson

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/MatrixScoringScheme.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>

#include <core/chemical/AA.hh>
#include <basic/database/open.hh>

#include <iostream>
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

static basic::Tracer tr( "core.sequence.MatrixScoringScheme" );

void MatrixScoringScheme::read_data( utility::io::izstream & input ) {
	std::string line;

	using utility::vector1;
	using namespace core::chemical;
	vector1< AA > order;
	while ( getline( input, line ) ) {
		if ( line.substr(0,1) == "#" ) continue; // skip comments

		std::istringstream line_stream( line );
		if ( line.substr(0,1) == " " ) { // header line
			char aa;
			while ( line_stream >> aa, !line_stream.fail() ) {
				if ( oneletter_code_specifies_aa( aa ) ) {
					order.push_back( aa_from_oneletter_code(aa) );
				} else {
					order.push_back( (core::chemical::AA) 0 ); //"Invalid" sentinel, to keep spacing correct
				}
			}
			vector1< Real > dummy( (core::Size) num_aa_types, 0.0 );
			scoring_matrix_.resize( (core::Size) num_aa_types, dummy );
		} else {
			char aa_name;
			line_stream >> aa_name;
			if ( !oneletter_code_specifies_aa( aa_name ) ) continue; // skip non-AA lines

			AA aa = aa_from_oneletter_code(aa_name);
			for ( vector1< AA >::const_iterator it = order.begin(),
					end = order.end(); it != end; ++it
					) {
				Real score;
				line_stream >> score;

				if ( line_stream.fail() ) {
					std::string message = "Error reading line " + line + '\n';
					utility_exit_with_message( message );
				}

				if ( (Size) aa <= order.size() && *it != 0 && (Size) *it <= order.size() ) {
					scoring_matrix_[ aa ][ *it ] = score;
				}

			}
		}
	} // while( getline( input, line ) )
} // read_data

/// @brief Read an alignment matrix from the given filename using the NCBI BLOSUM format
/// for matrices.
void MatrixScoringScheme::read_from_file( utility::file::FileName const & fn ) {
	utility::io::izstream input( fn );
	if ( !input ) {
		utility_exit_with_message(
			"ERROR: Unable to open MatrixScoringScheme file!" + std::string(fn)
		);
	}
	read_data( input );
} // read_from_file

/// @brief Read an alignment matrix from the given database filename using the
/// NCBI BLOSUM format for matrices.
void MatrixScoringScheme::read_from_database( std::string name ) {
	read_from_file( basic::database::full_name( "sequence/substitution_matrix/" + name ) );
}

/// @brief Get the values for amino acid aa, in Rosetta aa order.
utility::vector1< Real > MatrixScoringScheme::values_for_aa( core::chemical::AA aa ) {
	return scoring_matrix_[ aa ];
}

/// @brief Get the values for amino acid aa, in Rosetta aa order.
utility::vector1< Real > MatrixScoringScheme::values_for_aa( char aa ) {
	if ( core::chemical::oneletter_code_specifies_aa(aa) ) {
		return scoring_matrix_[ core::chemical::aa_from_oneletter_code(aa) ];
	} else {
		utility::vector1< Real > retval( scoring_matrix_[core::chemical::aa_ala].size() , 0.0 );
		return retval;
	}
}

utility::vector1< utility::vector1< Real > > MatrixScoringScheme::scoring_matrix() const {
	return scoring_matrix_;
}

Real MatrixScoringScheme::score(
	SequenceOP seq1,
	SequenceOP seq2,
	Size pos1,
	Size pos2
) {
	runtime_assert( pos1 <= seq1->length() );
	runtime_assert( pos2 <= seq2->length() );

	core::chemical::AA
		aa1( core::chemical::aa_from_oneletter_code( (*seq1)[pos1] ) ),
		aa2( core::chemical::aa_from_oneletter_code( (*seq2)[pos2] ) );

	if ( aa1 == core::chemical::aa_unk || aa2 == core::chemical::aa_unk ) {
		// likely a non-canonical aa in sequence
		tr.Error  << "returning score of zero for comparing amino acids "
			<< (*seq1)[pos1] << " and " << (*seq2)[pos2] << std::endl;
		return 0;
	}

	return 0.1 * scoring_matrix_[ aa1 ][ aa2 ];
}

} // sequence
} // core

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MatrixScoringScheme.hh
/// @brief class definition for a given scoring scheme for an alignment.
/// @details Simply based on comparing single characters from two protein
/// sequences, along with affine gap penalties of the form penalty = A + Bk, where
/// A represents the penalty for starting a gap, and B represents the penalty for
/// extending a previously opened gap by k characters.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_MatrixScoringScheme_hh
#define INCLUDED_core_sequence_MatrixScoringScheme_hh

#include <core/types.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/chemical/AA.hh>

#include <utility/io/izstream.fwd.hh>
#include <utility/file/FileName.fwd.hh>

#include <utility/vector1_bool.hh>


namespace core {
namespace sequence {

class MatrixScoringScheme : public ScoringScheme {

public:

	MatrixScoringScheme() {
		// gap open and gap extend set to default parameters
		// for NCBI BLASTP.
		gap_open  ( -11 );
		gap_extend( -1  );
		init_type();
	}

	/// @brief ctor
	MatrixScoringScheme(
		Real open,
		Real extend,
		utility::file::FileName const & fn
	)
	{
		read_from_file( fn );
		gap_open  ( open );
		gap_extend( extend );
		init_type();
	}

	MatrixScoringScheme(
		Real open,
		Real extend,
		utility::vector1< utility::vector1< core::Real > > matrix
	) : scoring_matrix_( matrix ) {
		gap_open  ( open );
		gap_extend( extend );
		init_type();
	}

	/// @brief returns owning pointer to a new object with a deep copy of
	/// this object's values.
	ScoringSchemeOP clone() const override {
		return ScoringSchemeOP( new MatrixScoringScheme(
			gap_open(),
			gap_extend(),
			scoring_matrix()
			) );
	}

	/// @brief dtor
	~MatrixScoringScheme() override = default;

	/// @brief Read an alignment matrix from the given database filename using the
	/// NCBI BLOSUM format for matrices.
	void read_from_database( std::string name="BLOSUM62" );
	/// @brief Read an alignment matrix from the given filename using the
	/// NCBI BLOSUM format for matrices.
	void read_from_file( utility::file::FileName const & fn ) override;
	void read_data( utility::io::izstream & input ) override;

	/// @brief Get the values for amino acid aa, in Rosetta aa order.
	utility::vector1< Real > values_for_aa( char aa );
	/// @brief Get the values for amino acid aa, in Rosetta aa order.
	utility::vector1< Real > values_for_aa( core::chemical::AA aa );
	utility::vector1< utility::vector1< Real > > scoring_matrix() const;

	Real score( SequenceOP seq1, SequenceOP seq2, Size pos1, Size pos2 ) override;

private:
	void init_type() {
		type( "Matrix" );
	}

	utility::vector1< utility::vector1< Real > > scoring_matrix_;
}; // class MatrixScoringScheme

} // sequence
} // core

#endif

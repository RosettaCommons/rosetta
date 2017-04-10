// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    AAMatrix.hh

/// @brief   Declarations and simple accessor/mutator definitions for AAMatrix.
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_AAMatrix_HH
#define INCLUDED_protocols_mean_field_AAMatrix_HH

// Unit header
#include <protocols/mean_field/AAMatrix.fwd.hh>

// Package headers
#include <protocols/mean_field/RotMatrix.fwd.hh>
#include <protocols/mean_field/AAProb.hh>
// Project headers
#include <core/chemical/AA.hh>


// Utility headers
#include <utility/vector1.hh>
#include <protocols/mean_field/jagged_array.hh>
#include <core/types.hh>

// Numeric headers


// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {
/// @details
/// @remarks
class AAMatrix : public protocols::mean_field::jagged_array< AAProb > {
public:
	// Standard methods ////////////////////////////////////////////////////////
	/// @brief Standard constructor
	AAMatrix();

	/// @brief Constructor to build AAMatrix from a file
	AAMatrix( std::istream & aa_matrix_file );

	/// @brief Constructs a AAMatrix from a RotMatrix (standard way of constructing AAMatrix)
	AAMatrix( RotMatrix const & rm,
		protocols::mean_field::jagged_array< core::Real > em, core::Real temp );

	/// @brief  Copy constructor
	AAMatrix( AAMatrix const & object_to_copy );

	/// @brief Assignment operator
	AAMatrix & operator=( AAMatrix const & object_to_copy );

	/// @brief Destructor
	~AAMatrix();

	/// @brief calculates vector of cosine distances between vectors of this AAMatrix and a second AAMatrix
	/// @details number of elements is the number of positions + 1.  last element is the cosine distance between entire matrices
	utility::vector1< core::Real > cosine_distance ( AAMatrix const & aa_matrix ) const;

	/// @brief calculates vector of frob distances between vectors of this AAMatrix and a second AAMatrix
	/// @details number of elements is the number of positions + 1.  last element is the frob dist between entire matrices
	utility::vector1< core::Real > frob_distance ( AAMatrix const & aa_matrix ) const;

	/// @brief calculates vector of average absolute differences between vectors of this AAMatrix and a second AAMatrix
	/// @details number of elements is the number of positions + 1.  last element is the average abs diff between entire matrices
	utility::vector1< core::Real > ave_abs_diff ( AAMatrix const & aa_matrix ) const;

	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief  Generate string representation of AAMatrix for debugging purposes.
	void show( std::ostream & output=std::cout ) const;

	/// @brief outputs AA Matrix in transfac format for dumping purposes
	void dump_transfac( std::string const & filename )  const;

	/// @brief returns true if all probabilities are 0
	bool empty() const;

	/// @brief Insertion operator (overloaded so that AAMatrix can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, AAMatrix const & object_to_output );

	// Accessors/Mutators
	/// @brief returns vector of current amino acids at each position (if current rotamer is set for that position)
	inline
	utility::vector1< core::Size > const & curr_aa() const
	{
		return curr_aa_;
	}

	/// @brief returns vector of probabilities of current amino acids at each position (if current rotamer is set for that position)
	utility::vector1< AAProb > curr_aa_prob() const;

	/// @brief builds AAMatrix from RotMatrix rm using private method init
	inline
	void build_aa_matrix( RotMatrix const & rm,
		protocols::mean_field::jagged_array< core::Real > em,
		core::Real temp )
	{
		init ( rm, em, temp );
	}

	/// @brief builds AAMatrix from transfac file using private method init
	inline
	void build_aa_matrix( std::istream & aa_matrix_file )
	{
		init ( aa_matrix_file );
	}

private:
	// Private methods /////////////////////////////////////////////////////////
	/// @brief Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data( AAMatrix & object_to_copy_to, AAMatrix const & object_to_copy_from );

	/// @brief private method called by constructors to initialize the AAMatrix from a RotMatrix
	void init( protocols::mean_field::RotMatrix const & rm,
		protocols::mean_field::jagged_array< core::Real > em, core::Real temp );

	/// @brief private method called by constructors to initialize the AAMatrix from a transfac file
	void init( std::istream & aa_matrix_file );

	/// @brief parse line of transfac file to add a vector of probabilities for one position to the AAMatrix
	utility::vector1< AAProb > parse_aa_matrix_line( utility::vector1< std::string > const & tokens,
		utility::vector1< core::chemical::AA > const & aa_names );

	/// @brief parse header (AA codes) line of the AAMatrix
	void parse_aa_line( utility::vector1< std::string > const & tokens, utility::vector1< core::chemical::AA > & aa_names );

	// Private data ////////////////////////////////////////////////////////////
	utility::vector1 < core::Size > curr_aa_;
};  // class AAMatrix

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_AAMatrix_HH

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    RotMatrix.hh

/// @brief   Declarations and simple accessor/mutator definitions for RotMatrix (a conformational matrix used for mean-field theory calculations).
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_RotMatrix_HH
#define INCLUDED_protocols_mean_field_RotMatrix_HH

// Unit header
#include <protocols/mean_field/RotMatrix.fwd.hh>

// Package headers
#include <protocols/mean_field/RotProb.hh>

// Project headers
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <protocols/mean_field/jagged_array.hh>
#include <core/types.hh>

// Numeric headers


// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {

/// @details derived from a jagged_array of RotProbs
/// @remarks used in MeanField family of classes
class RotMatrix : public protocols::mean_field::jagged_array< RotProb > {
public:
	// Standard methods ////////////////////////////////////////////////////////

	/// @brief default constructor, can be useful
	RotMatrix();

	/// @brief Standard constructor
	RotMatrix( core::Size const option, core::pack::rotamer_set::RotamerSetsOP rs );

	/// @brief  Copy constructor
	RotMatrix( RotMatrix const & object_to_copy );

	/// @brief Assignment operator
	RotMatrix & operator=( RotMatrix const & object_to_copy );

	/// @brief Destructor
	~RotMatrix();


	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief  Generate string representation of RotMatrix for debugging purposes.
	void show( std::ostream & output=std::cout ) const;

	/// @brief Insertion operator (overloaded so that RotMatrix can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, RotMatrix const & object_to_output );


	// Accessors/Mutators

	/// @brief returns a vector of the current rotamers at each position
	inline
	utility::vector1< core::Size > const & curr_rot() const
	{
		return curr_rot_;
	}

	/// @brief sets the curr_rot vector
	inline
	void curr_rot( utility::vector1< core::Size > cr )
	{
		curr_rot_ = cr;
	}

	/// @brief returns a vector of bools designating whether each position allows design or not
	inline
	utility::vector1< bool > const & is_designed() const
	{
		return is_designed_;
	}

	/// @brief sets a vector of bools designating whether each position allows design or not
	inline
	void is_designed( utility::vector1< bool > id )
	{
		is_designed_ = id;
	}

	/// @brief returns a bool designating whether a specified position allows design
	/// @param [in] pos - position to be examined within the RotMatrix
	/// @remarks RotMatrix numbering may or may not match pose numbering
	inline
	bool is_designed( core::Size pos ) const
	{
		return is_designed_[ pos ];
	}

	/// @brief returns the number of positions that allow design
	core::Size n_designed() const;

	/// @brief returns a vector of RotProbs corresponding to the current rotamer at each position
	utility::vector1 < RotProb > curr_rot_prob() const;

	/// @brief calls init so as to reinitialize the RotMatrix with the given params
	void build_rot_matrix( core::Size const option, core::pack::rotamer_set::RotamerSetsOP rs );

private:
	// Private methods /////////////////////////////////////////////////////////
	/// @brief Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data( RotMatrix & object_to_copy_to, RotMatrix const & object_to_copy_from );

	/// @brief private method called by constructors to initialize the RotMatrix
	void init( core::Size const option, core::pack::rotamer_set::RotamerSetsOP rs );

	// Private data ////////////////////////////////////////////////////////////
	utility::vector1 < core::Size > curr_rot_; //length of positions in RotMatrix
	utility::vector1 < bool > is_designed_; //length of positions in RotMatrix
};  // class RotMatrix
}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_RotMatrix_HH

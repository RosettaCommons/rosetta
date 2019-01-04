// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//////////////////////////////////////////////////////////////////////
///
/// @brief Generic base class for the MathNTensor class.  Since the MathNTensor class takes
/// a type AND a dimensionality as template arguments, it's not possible to have a generic
/// pointer to a MathNTensor of arbitrary dimensionality.  The base class allows this.
///
/// @author Vikram K. Mulligan (vmullig@uw.edu).
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_MathNTensorBase_hh
#define INCLUDED_numeric_MathNTensorBase_hh

// Package headers
#include <numeric/types.hh>
#include <numeric/MathNTensorBase.fwd.hh>
#include <numeric/MathNTensor.fwd.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// C++ headers
#include <math.h>
#include <iostream>
#include <memory>

namespace numeric {

template< class T >
class MathNTensorBase
{
public:
	typedef numeric::Size Size;

	/// @brief Default constructor.
	///
	MathNTensorBase() :
		dimensionality_(0)
	{}

	/// @brief Virtual destructor needed for polymorphism.
	///
	virtual ~MathNTensorBase() {}

	/// @brief Pure virtual clone function.
	/// @details Creates a copy of this object and returns an owning pointer to the copy.  Must be
	/// implemented by derived classes.
	virtual MathNTensorBaseOP< T > clone() const = 0;

	/// @brief Costructor with dimensionality value.
	///
	MathNTensorBase( Size const dimensionality_in ) :
		dimensionality_(dimensionality_in)
	{}

	/// @brief Get the dimensionality of derived classes.
	/// @details Will need to store this in order to cast pointers to the appropriate type for the derived class.
	Size dimensionality() const { return dimensionality_; }

	/// @brief Is the given coordinate a valid position in the tensor?
	/// @details Returns false if out of range.  Pure virtual; must be implemented by derived classes.
	virtual bool is_valid_position( utility::vector1< Size > const & position ) const = 0;

	/// @brief Set a value in a tensor.
	/// @details Note that bounds-checking only occurs in debug builds!
	/// @note Pure virtual.  Must be implemented by derived classes.
	virtual void set_value( utility::vector1< Size > const &position, T const &value_in ) = 0;

	/// @brief Get a value from a tensor.
	/// @details Note that bounds-checking only occurs in debug builds!
	/// @note Pure virtual.  Must be implemented by derived classes.
	virtual T const & value( utility::vector1< Size > const &position ) const = 0;

	/// @brief Get the number of bins for the nth dimension.
	/// @details Pure virtual.  Must be implemented by derived classes.
	virtual Size n_bins( Size const dimension ) const = 0;

	/// @brief Get the minimum value stored in this tensor.
	/// @details pure virtual.  Must be implemented by derived classes.
	virtual T min() const = 0;

	/// @brief Get the maximum value stored in this tensor.
	/// @details pure virtual.  Must be implemented by derived classes.
	virtual T max() const = 0;

protected:

	/// @brief Lets the derived class set the dimensionality stored in the base class.
	/// @details Will need to store this in order to cast pointers to the appropriate type for the derived class.
	void set_dimensionality( Size const dimensionality_in ) { dimensionality_ = dimensionality_in; }

private:
	/// @brief The dimensionality of the MathNTensor/
	///
	Size dimensionality_ = 0;

};

template< class T >
MathNTensorBaseOP< T >
deep_copy( MathNTensorBase< T > const & source ) {
	return source.clone();
}

/// @brief Explicit template instantiation, apparently needed for PyRosetta
template MathNTensorBaseOP< Real > deep_copy( MathNTensorBase< Real > const &  );

}//end namespace numeric


#endif

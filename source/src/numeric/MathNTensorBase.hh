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

	/// @brief Costructor with dimensionality value.
	///
	MathNTensorBase( Size const dimensionality_in ) :
		dimensionality_(dimensionality_in)
	{}

	/// @brief Destructor.
	///
	~MathNTensorBase() {}

	/// @brief Get the dimensionality of derived classes.
	/// @details Will need to store this in order to cast pointers to the appropriate type for the derived class.
	Size dimensionality() const { return dimensionality_; }

protected:

	/// @brief Lets the derived class set the dimensionality stored in the base class.
	/// @details Will need to store this in order to cast pointers to the appropriate type for the derived class.
	void set_dimensionality( Size const dimensionality_in ) { dimensionality_ = dimensionality_in; }

private:
	/// @brief The dimensionality of the MathNTensor/
	///
	Size dimensionality_=0;

};


}//end namespace numeric


#endif

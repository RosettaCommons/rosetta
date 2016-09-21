// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MathNTensor.fwd.hh
/// @brief  Owning pointer declarations for MathNTensor template class.
/// @details Note that if you're using an owning pointer to a MathNTensor, you need to either know at pointer declaration time
/// the dimensionality of the MathNTensor, or you need to be using an owning pointer to the base class (MathNTensorBase).
/// @author Vikram K. Mulligan (vmullig@uw.edu).


#ifndef INCLUDED_numeric_MathNTensor_fwd_hh
#define INCLUDED_numeric_MathNTensor_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <numeric/types.hh>

namespace numeric {

// Forward declaration
template< class T, numeric::Size N > class MathNTensor;

// Owning pointer
template< class T, numeric::Size N >
using MathNTensorOP = utility::pointer::shared_ptr< MathNTensor< T,N > >; //Vikram is reluctantly using a C++11 feature.  Mrph.

// Const-access owning pointer
template< class T, numeric::Size N >
using MathNTensorCOP = utility::pointer::shared_ptr< MathNTensor< T,N > const >;

} // namespace numeric

#endif // INCLUDED_numeric_MathNTensor_fwd_hh

